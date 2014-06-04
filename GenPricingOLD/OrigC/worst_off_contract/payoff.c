/* Monte-Carlo pricing code generated by LexiFi's MLFi compiler
   Compiler version: 3.12.1+dev4 (2010-09-03).
   Do not edit by hand! */
#include <mlfi_pricing.h>
static void trajectory(const mlfi_nmc_model* model) {
  double **underlyings = model->underlyings;
  double x50;

  if ((1. <= fmin((underlyings[0][1] / 11840.), fmin((underlyings[0][2] / 1200.), (underlyings[0][0] / 3758.05))))) goto L21;
  if ((1. <= fmin((underlyings[1][1] / 11840.), fmin((underlyings[1][2] / 1200.), (underlyings[1][0] / 3758.05))))) goto L22;
  if ((1. <= fmin((underlyings[2][1] / 11840.), fmin((underlyings[2][2] / 1200.), (underlyings[2][0] / 3758.05))))) goto L23;
  if ((1. <= fmin((underlyings[3][1] / 11840.), fmin((underlyings[3][2] / 1200.), (underlyings[3][0] / 3758.05))))) goto L24;
  x50=fmin((underlyings[4][1] / 11840.), fmin((underlyings[4][2] / 1200.), (underlyings[4][0] / 3758.05)));
  model->notify_cash_flow(model, 0, 1000., 4 /*2017-02-03, 2017-01-27, EUR*/);
  if ((1. <= x50)) goto L25;
  if ((0.75 < x50)) return;
  model->notify_cash_flow(model, 0, (-(1000. * (1. - x50))), 4 /*2017-02-03, 2017-01-27, EUR*/); return;
L25:
  model->notify_cash_flow(model, 0, 750., 4 /*2017-02-03, 2017-01-27, EUR*/); return;
L24:
  model->notify_cash_flow(model, 0, 1600., 3 /*2016-02-03, 2016-01-27, EUR*/); return;
L23:
  model->notify_cash_flow(model, 0, 1450., 2 /*2015-02-03, 2015-01-27, EUR*/); return;
L22:
  model->notify_cash_flow(model, 0, 1300., 1 /*2014-02-03, 2014-01-27, EUR*/); return;
L21:
  model->notify_cash_flow(model, 0, 1150., 0 /*2013-02-01, 2013-01-27, EUR*/);
  return;

}
static void init(const mlfi_nmc_model* model) {
}
static mlfi_date state_dates[5]=
{17398800 /*2013-01-27*/, 17924400 /*2014-01-27*/, 18450000 /*2015-01-27*/,
 18975600 /*2016-01-27*/, 19502640 /*2017-01-27*/};
static mlfi_cash_flow cash_flows[5]=
{{17406000 /*2013-02-01*/, 17398800 /*2013-01-27*/, EUR},
 {17934480 /*2014-02-03*/, 17924400 /*2014-01-27*/, EUR},
 {18460080 /*2015-02-03*/, 18450000 /*2015-01-27*/, EUR},
 {18985680 /*2016-02-03*/, 18975600 /*2016-01-27*/, EUR},
 {19512720 /*2017-02-03*/, 19502640 /*2017-01-27*/, EUR}};
static char* std_underlyings[3]=
{"DJ_Eurostoxx_50", "Nikkei_225", "SP_500"};
static mlfi_nmc_pricing_code mc_info ={
 /*nb_contracts*/1,
 5, state_dates,
 5, cash_flows,
 3, std_underlyings,
 0, /*notify_parameters*/NULL,
 0, /*custom_pricers*/NULL,
 0, /*deterministic_pricers*/NULL,
 /*only_initial_date_pricing*/0,
 /*size_of_pricing_data*/0,
 /*max_valuation_date*/168307199 /*2299-12-31T23:59:00*/,
 &init,
 &trajectory
};
MLFI_EXPORT mlfi_nmc_pricing_code* get_pricing_code(){return &mc_info;}
