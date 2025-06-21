/*
 * ETAS catalog generator include.
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

#ifndef ETASCATGEN_ETASCATGEN_HPP
#define ETASCATGEN_ETASCATGEN_HPP

#include <cyantities/unit.hpp>
#include <cyantities/quantitywrap.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/frequency.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

namespace etascatgen {

namespace bu = boost::units;


using Time = bu::quantity<bu::si::time, double>;
using Frequency = bu::quantity<bu::si::frequency, double>;
using Scalar = bu::quantity<bu::si::dimensionless, double>;

/*
 * Earthquake with magnitude and occurrence time (no spatial information):
 */

void ETAS_generate_catalog_M_t(
    const cyantities::QuantityWrapper& mu_0,
    double Mmin,
    double Mmax,
    double beta,
    double alpha,
    double p,
    const cyantities::QuantityWrapper& c,
    double offspring_fraction,
    const size_t N_skip,
    size_t seed,
    cyantities::QuantityWrapper& Mi,
    cyantities::QuantityWrapper& ti
);

}

#endif