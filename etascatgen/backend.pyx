# ETAS catalog generator backend wrapper.
#
# Author: Malte J. Ziebarth (mjz.science@fmvkb.de)
#
# Copyright (C) 2025 Malte J. Ziebarth
#
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
# the European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Licence is distributed on an "AS IS" basis,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and
# limitations under the Licence.

from cyantities.quantity cimport Quantity, QuantityWrapper

cdef extern from "etascatgen/etascatgen.hpp" namespace "etascatgen" nogil:
    void ETAS_generate_catalog_M_t(
        const QuantityWrapper& mu_0,
        double Mmin,
        double Mmax,
        double beta,
        double alpha,
        double p,
        const QuantityWrapper& c,
        double offspring_fraction,
        size_t N_skip,
        size_t seed,
        QuantityWrapper& Mi,
        QuantityWrapper& ti
    ) except+





def generate_catalog_M_t(
        size_t N,
        Quantity mu_0,
        double Mmin,
        double Mmax,
        double beta,
        double alpha,
        double p,
        Quantity c,
        double offspring_fraction,
        size_t N_skip,
        size_t seed = 198372
    ):
    assert mu_0._is_scalar

    cdef Quantity Mi = Quantity.zeros(N, '1')
    cdef Quantity ti = Quantity.zeros(N, 's')

    ETAS_generate_catalog_M_t(
        mu_0.wrapper(),
        Mmin,
        Mmax,
        beta,
        alpha,
        p,
        c.wrapper(),
        offspring_fraction,
        N_skip,
        seed,
        Mi.wrapper(),
        ti.wrapper()
    )

    return Mi, ti