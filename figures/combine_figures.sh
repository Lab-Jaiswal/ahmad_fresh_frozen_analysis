#!/usr/bin/env bash

set -o xtrace -o nounset -o pipefail -o errexit

qpdf --empty --pages ../figure1/figure1.pdf \
    ../figure2/figure2.pdf 1-z \
    ../figure3/figure3.pdf 1-z \
    ../figure4/figure4.pdf 1-z \
    ../figure_s1/figure_s1.pdf 1-z \
    ../figure_s3/figure_s3.pdf 1-z \
    ../figure_s4/figure_s4.pdf 1-z \
    ../figure_s5/figure_s5.pdf 1-z \
    ../figure_s6/figure_s6.pdf 1-z \
    ../figure_s7/figure_s7.pdf 1-z \
    -- figures.pdf

qpdf --empty --pages ../figure1/figure1.pdf \
    ../figure2/figure2.pdf 1-z \
    ../figure3/figure3.pdf 1-z \
    ../figure4/figure4.pdf 1-z \
    ../figure_s1/figure_s1_caption.pdf 1-z \
    ../figure_s3/figure_s3_caption.pdf 1-z \
    ../figure_s4/figure_s4_caption.pdf 1-z \
    ../figure_s5/figure_s5_caption.pdf 1-z \
    ../figure_s6/figure_s6_caption.pdf 1-z \
    ../figure_s7/figure_s7_caption.pdf 1-z \
    -- figures_supp_captions.pdf
