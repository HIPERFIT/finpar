import pylab as P

"""
swaption: maturity, swap_term, frequency (months)
'quotes': Black '76 (implied) volatility
"""

swaption_quotes = [
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  1, 'swap_frequency': 6}, 1.052),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.81485),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.6165),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.46995),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.38295),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.3325),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.3016),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.2815),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.26435),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.2496),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.2516),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.28835),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.27155),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  1, 'swap_frequency': 6}, 0.23465),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.61445),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.54805),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.46795),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.3919),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.3434),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.3083),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.28655),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.2697),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.25775),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.2443),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.26495),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.28195),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.26845),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  2, 'swap_frequency': 6}, 0.20995),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.5835),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.49255),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.42825),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.3695),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.329),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.3022),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.28165),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.26615),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.25485),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.24375),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.2718),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.28135),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.26865),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  3, 'swap_frequency': 6}, 0.2131),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.5415),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.46235),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.403),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.3559),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.3232),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.29675),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.27715),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.26385),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.254),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.2454),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.27845),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.2821),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.2678),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  4, 'swap_frequency': 6}, 0.2131),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.517),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.446),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.3903),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.34755),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.3166),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.29305),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.2745),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.2639),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.2534),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.2499),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.28315),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.2825),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.277),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  5, 'swap_frequency': 6}, 0.21175),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.478),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.42105),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.37715),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.3378),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.311),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.2895),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.2745),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.264),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.2573),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.25475),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.28815),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.28195),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.26015),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  6, 'swap_frequency': 6}, 0.2097),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.452),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.4074),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.368),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.3307),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.30645),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.2877),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.27475),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.2664),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.26155),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.26035),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.292),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.2825),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.25685),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  7, 'swap_frequency': 6}, 0.2081),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.43395),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.39445),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.35885),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.3281),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.30395),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.28745),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.2767),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.27065),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.26625),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.26625),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.2921),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.2814),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.25265),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  8, 'swap_frequency': 6}, 0.2083),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.42285),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.3857),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.3521),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.3239),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.30285),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.2895),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.2799),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.27485),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.2712),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.27205),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.29205),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.27855),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.24945),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year':  9, 'swap_frequency': 6}, 0.219),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.41765),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.38095),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.34795),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.3217),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.30365),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.2916),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.2842),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.27985),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.2769),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.2775),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.306),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.2763),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.2458),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year': 10, 'swap_frequency': 6}, 0.22),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.37905),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.35465),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.33505),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.31725),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.3008),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.29075),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.28365),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.2787),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.27385),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.2709),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.2689),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.24225),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.2096),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year': 15, 'swap_frequency': 6}, 0.18285),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.37975),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.3605),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.3407),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.321),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.3063),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.29315),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.28395),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.2777),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.27205),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.26675),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.24875),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.21735),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.1939),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year': 20, 'swap_frequency': 6}, 0.17205),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.38115),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.3627),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.34425),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.3222),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.3084),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.2941),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.28285),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.2751),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.2663),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.26055),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.2338),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.20735),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.1823),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year': 25, 'swap_frequency': 6}, 0.1686),
  ({'swaption_maturity_in_year':  1, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.38285),
  ({'swaption_maturity_in_year':  2, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.3633),
  ({'swaption_maturity_in_year':  3, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.34125),
  ({'swaption_maturity_in_year':  4, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.3188),
  ({'swaption_maturity_in_year':  5, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.30305),
  ({'swaption_maturity_in_year':  6, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.2888),
  ({'swaption_maturity_in_year':  7, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.2748),
  ({'swaption_maturity_in_year':  8, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.26725),
  ({'swaption_maturity_in_year':  9, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.25985),
  ({'swaption_maturity_in_year': 10, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.25165),
  ({'swaption_maturity_in_year': 15, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.2267),
  ({'swaption_maturity_in_year': 20, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.1989),
  ({'swaption_maturity_in_year': 25, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.18115),
  ({'swaption_maturity_in_year': 30, 'swap_term_in_year': 30, 'swap_frequency': 6}, 0.16355),
  ]

def maturities(swaption_quotes=swaption_quotes):
  l=list( set( map(lambda (swp,vol): swp['swaption_maturity_in_year'],swaption_quotes) ) )
  l.sort()
  return l

def swap_terms(swaption_quotes=swaption_quotes):
  l=list( set( map(lambda (swp,vol): swp['swap_term_in_year'],swaption_quotes) ) )
  l.sort()
  return l

def plot3Dwireframe(x,y,z,ax_x,ax_y,xlabel=None,ylabel=None,zlabel=None,title=None,fig=None):
  import pylab as P
  import numpy as N
  from mpl_toolkits.mplot3d import Axes3D, axes3d
  #
  if fig == None: fig=P.figure()
  ax = fig.add_subplot(111, projection='3d')
  _x,_y=N.meshgrid(ax_x,ax_y)
  ax.plot_wireframe(_x,_y,N.array(z).reshape((len(ax_x),len(ax_y))))
  ax.plot(x,y,z,'bo')
  #
  if xlabel: ax.set_xlabel(xlabel)
  if ylabel: ax.set_ylabel(ylabel)
  if zlabel: ax.set_zlabel(zlabel)
  if title: ax.set_title(title)
  #
  return ax

def swaption_plot(swaption_quotes=swaption_quotes):
  import numpy as N
  #
  xlabel='swaption maturity in year'
  ylabel='swap term in year'
  zlabel='swaption quote'
  title="quotes: implied volatilities"
  quotes=N.array( map( lambda swp: (swp[0][xlabel.replace(" ","_")],swp[0][ylabel.replace(" ","_")],swp[1]) , swaption_quotes) )
  #
  _m=maturities(swaption_quotes=swaption_quotes)
  _t=swap_terms(swaption_quotes=swaption_quotes)
  plot3Dwireframe(x=quotes[:,0],y=quotes[:,1],z=quotes[:,2],ax_x=_m,ax_y=_t,xlabel=xlabel,ylabel=ylabel,zlabel=zlabel,title=title)

if __name__=='__main__':
  swaption_plot(swaption_quotes=swaption_quotes)
  P.show()
