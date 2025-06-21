/*
 * ETAS catalog generator without spatial information.
 *
 * Author: Malte J. Ziebarth (mjz.science@fmvkb.de)
 *
 * Copyright (C) 2025 Malte J. Ziebarth
 *
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
 * the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 */

#include <etascatgen/etascatgen.hpp>
#include <ranges>
#include <numeric>
#include <random>
#include <optional>
#include <queue>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/tools/roots.hpp>

#include <iostream>


namespace etascatgen {

/*
 * Parameters used in this implementation
 * ======================================
 *
 * FK : Frequency `K / Tref^p` derived from K of Ogata (1988)
 *      and a reference time scale Tref.
 *      We use this instead of K itself to avoid fractional
 *      units.
 */

struct Process_M_t {
    Frequency mu_0;
    Time Tref;
    Time c;
    double beta;
    double Mr;
    double p;
    double ln_p;
    Frequency FK;

    Process_M_t(
        Frequency mu_0,
        Time Tref,
        Time c,
        double beta,
        double Mr,
        double p,
        double Mmin,
        double Mmax,
        double offspring_fraction
    ) : mu_0(mu_0), Tref(Tref), c(c), beta(beta),
        Mr(Mr), p(p), ln_p(std::log(p)),
        FK(
            critical_FK(Mmin, Mmax, p, c, Tref, beta, Mr)
            * offspring_fraction
        )
    {}

private:
    static Frequency critical_FK(
        double Mmin,
        double Mmax,
        double p,
        Time c,
        Time Tref,
        double beta,
        double Mr
    )
    {
        /*
         * A translation of the Python code:
         *
         * (p-1) * c ** (p-1)
         *   * (1 - np.exp(-beta * (Mmax - Mmin)))
         *   / (beta * np.exp(beta * (Mmin - Mr))*(Mmax - Mmin))
         *   / Tref ** p
         */
        return (p - 1.0) *
            std::exp(p * std::log(c / Tref)) / c
            * (1.0 - std::exp(-beta * (Mmax - Mmin)))
            / (beta * std::exp(beta * (Mmin - Mr)) * (Mmax - Mmin));
    }
};


// static Frequency single_rate(
//     Time t,
//     Time ti,
//     Scalar Mi,
//     const Process_M_t& process
// )
// {
//     return process.FK * std::exp(
//         -process.beta * (Mi - process.Mr)
//         + process.ln_p * std::log(
//             process.Tref / (t - ti + process.c)
//         )
//     );
// }


/*
 * Compute the time of the next descendant
 */
static double f(double M, const Process_M_t& process)
{
    return std::exp(process.beta * (M - process.Mr));
}

// static double Lambda_i(
//     Time ti,
//     Time tl,
//     Time tr,
//     double Mi,
//     const Process_M_t& process
// )
// {
//     /*
//      * FK = K / Tref ** p
//      * Thereby:
//      *    Tref * FK * ((tr - ti + c)/Tref) ** (1-p)
//      *    = Tref * (K * Tref ** -p) * (tr - ti + c) ** (1-p)
//      *      * Tref ** (p - 1)
//      *    = K * (tr - ti + c) ** (1-p)
//      */
//     double _1mp = 1.0 - process.p;
//     return f(Mi, process) * process.Tref * process.FK / _1mp * (
//         std::pow((tr - ti + process.c) / process.Tref, _1mp)
//         - std::pow((tl - ti + process.c) / process.Tref, _1mp)
//     );
// }

static double Lambda_i_oo(
    Time ti,
    Time tl,
    double Mi,
    const Process_M_t& process
)
{
    /*
     * FK = K / Tref ** p
     * Thereby:
     *    Tref * FK * ((tr - ti + c)/Tref) ** (1-p)
     *    = Tref * (K * Tref ** -p) * (tr - ti + c) ** (1-p)
     *      * Tref ** (p - 1)
     *    = K * (tr - ti + c) ** (1-p)
     */
    double _1mp = 1.0 - process.p;
    return -f(Mi, process) * process.Tref * process.FK / _1mp *
        std::pow((tl - ti + process.c) / process.Tref, _1mp);
}


static std::optional<Time> next_single_occurrence(
    double q,
    Time ti,
    double Mi,
    Time tl,
    const Process_M_t& process
)
{
    /* Early exit if no occurrence in finite time: */
    if (q <= std::exp(-Lambda_i_oo(ti, tl, Mi, process)))
        return std::optional<Time>();

    /*
     * Note here that we extract a factor Tref ** (1 - p)
     * from the outer logarithm.
     * The first summand has just the right exponent! So
     * we can simply divide by Tref.
     * The second summand is more tricky. However, note that
     *    K = FK * Tref ** p,
     * so that
     *    (1/K) / Tref ** (1 - p)
     *       = 1 / (FK * Tref ** p * Tref ** (1 - p))
     *       = 1 / (FK * Tref)
     */
    double _1mp = 1.0 - process.p;
    return ti - process.c + process.Tref * std::exp(
        1.0 / _1mp * std::log(
            std::pow((tl - ti + process.c) / process.Tref, _1mp)
            - _1mp / (f(Mi, process) * process.FK * process.Tref) * std::log(q)
       )
    );
}


static Time next_background_occurrence(
    double q,
    Time tl,
    const Process_M_t& process
)
{
    return tl - std::log(q) / process.mu_0;
}


static double draw_magnitude(
    double q,
    double Mmin,
    double Mmax,
    double beta
)
{
    return Mmin - std::log(
            1.0 - q * (1.0 - std::exp(-beta * (Mmax - Mmin)))
    ) / beta;
}


/*
 * This structure holds the components of the Hawkes process
 * intensity:
 */
struct excitement_t {
    Time ti;
    double M;
    Time tnext;

    /*
     * Ordering for a priority queue that provides the
     * next event:
     */
    bool operator>(const excitement_t& other) const
    {
        return tnext < other.tnext;
    }

    bool operator<(const excitement_t& other) const
    {
        return tnext > other.tnext;
    }
};




void ETAS_generate_catalog_M_t(
    const cyantities::QuantityWrapper& mu_0,
    double Mmin,
    double Mmax,
    double beta,
    double p,
    const cyantities::QuantityWrapper& c,
    double Mr,
    double offspring_fraction,
    const size_t N_skip,
    size_t seed,
    cyantities::QuantityWrapper& Mi,
    cyantities::QuantityWrapper& ti
)
{
    /* Sanity: */
    if (Mmin >= Mmax)
        throw std::runtime_error("Mmin >= Mmax");
    //if (beta < )
    if (p <= 1.0)
        throw std::runtime_error("p <= 1");
    if (offspring_fraction >= 1.0)
        throw std::runtime_error("Instable process (offspring ratio > 1)");
    else if (offspring_fraction < 0.0)
        throw std::runtime_error("Offspring ratio needs to be non-negative.");

    const size_t N = Mi.size();
    if (ti.size() != N)
        throw std::runtime_error("Size of M and t not compatible");


    /* Current number of earthquakes generated: */
    size_t n = 0;

    /* Normalization: */
    constexpr Time Tref = 1.0 * bu::si::seconds;

    Process_M_t process(
        mu_0.get<Frequency>(),
        Tref,
        c.get<Time>(),
        beta,
        Mr,
        p,
        Mmin,
        Mmax,
        offspring_fraction
    );

    /* Init the RNG: */
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);


    /* Time and magnitude of the current event: */
    Time t = 0.0 * bu::si::seconds;
    double M = std::numeric_limits<double>::quiet_NaN();

    /*
     * The next background occurrence:
     */
    Time next_bg = next_background_occurrence(
        uniform(rng),
        t,
        process
    );

    /*
     * A priority queue of future descendants of the intensity
     * components:
     */
    std::priority_queue<excitement_t> descendants;

    /*
     * The loop body:
     */
    auto next_event = [&]()
    {
        /*
         * Get the next occurrence time:
         */
        if (descendants.empty() || next_bg < descendants.top().tnext){
            /* Background event */
            t = next_bg;
            next_bg = next_background_occurrence(
                uniform(rng),
                t,
                process
            );
        } else {
            /* Descendant event. Pop it from the queue: */
            excitement_t event(descendants.top());
            descendants.pop();
            t = event.tnext;

            /* Check whether we generate a new descendant event from the
             * initial: */
            std::optional<Time> tnext(next_single_occurrence(
                uniform(rng),
                event.ti,
                event.M,
                t,
                process
            ));
            if (tnext){
                event.tnext = *tnext;
                descendants.push(event);
            }
        }


        /*
         * Get the next magnitude:
         */
        M = draw_magnitude(
            uniform(rng),
            Mmin,
            Mmax,
            beta
        );

        /*
         * Check whether this earthquake triggers another:
         */
        std::optional<Time> tnext(next_single_occurrence(
            uniform(rng),
            t,
            M,
            t,
            process
        ));
        if (tnext){
            descendants.push(
                excitement_t(
                    t,
                    M,
                    *tnext
                )
            );
        }
    };

    /*
     * First start up the state by throwing away the first couple of
     * earthquakes:
     */
    while (n < N_skip){
        next_event();
        ++n;
    }


    /*
     * Get some ranges and iterators for the working buffers:
     */
    auto M_out = Mi.iter<Scalar>();
    auto M_out_i = M_out.begin();
    auto t_out = ti.iter<Time>();
    auto t_out_i = t_out.begin();
    while (n < N+N_skip){
        /*
         * Generate next event:
         */
        next_event();

        /*
         * Save and advance the output iterators:
         */
        *t_out_i = t;
        *M_out_i = M;
        ++t_out_i;
        ++M_out_i;
        ++n;
    }
}


}