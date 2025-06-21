# ETASCatGen: An ETAS Earthquake Catalog Generator
The `etascatgen` Python module and C++ library generates space-agnostic earthquake
catalogs according to the Epidemic-Type Aftershock Sequence (ETAS) model of temporal
earthquake clustering. Specifically, ETASCatGen is based on the model described by
[Ogata (1988)](https://doi.org/10.2307/2288914). This model is a marked linear stationary
Hawkes process [(Hawkes, 1971)](doi.org/10.1111/j.2517-6161.1971.tb01530.x) with time-dependent intensity
```math
    \lambda(t) = \mu + \sum\limits_{t_i < t} f(M_i) g(t-t_i)\,,
```
where $t_i$ are the occurrence time of past earthquakes (the events), and $M_i$ are the magnitudes of the earthquakes (the markers).

In the epidemic-type model of aftershock excitement proposed by Ogata (1988), the rate of aftershock generation depends on the magnitude of generated earthquakes,
```math
    f(M) = \exp\big(\alpha(M - M_r)\big)\,,
```
where $M_r$ is the 'reference magnitude' (Ogata, 1988). In ETASCatGen, $M_r$ is absorbed in
the parameter describing the offspring fraction.

...

Once an earthquake occurs, the magnitude is drawn (independently) from a Gutenberg-Richter
distribution with an adjustable $b$-value confined to the magnitude interval
$M_\mathrm{min} \leq M \leq M_\mathrm{max}$,
```math
    \phi(M) = \frac{\beta \exp\big(-\beta (M - M_\mathrm{min})\big)}
                   {1 - \exp\big(-\beta(M_\mathrm{max} - M_\mathrm{min})\big)}
```
with $\beta = \ln(10) b$ (Ogata, 1988).

The ETASCatGen module takes a number of choices with regards to the parameters:
 - Without loss of generality, $M_r$ is chosen to be identical to $M_\mathrm{min}$.
 - The scale parameter $K$ is expressed as a fraction of $K^*$, the critical $K$
   at which the Hawkes process catastrophically self-resonates (that is, an expected
   offspring number of 1).