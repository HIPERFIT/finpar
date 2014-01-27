(***************************************************************************)
(*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               *)
(*                                                                         *)
(*  No part of this document may be reproduced or transmitted in any       *)
(*  form or for any purpose without the express permission of LexiFi SAS.  *)
(***************************************************************************)

(* $Id: g2pp.mfi 50183 2012-10-26 08:14:04Z vbegot $ *)

(** G2++ 2-factors affine interest rate model.  See Brigo and Mercurio
    (2006), chapter 4.2 p142 for more details. *)

open Date

type swap_frequency =
  | Freq_3M
  | Freq_6M
  | Freq_1Y

type swaption =
    {
     swaption_maturity_in_year: int;
     swap_frequency: swap_frequency;
     swap_term_in_year: int;
    }

type parameters =
    {
     g_a : float;
     g_b : float;
     g_sigma : float;
     g_nu : float;
     g_rho : float;
    }
(** A point in the calibration space of the G2++ model.
    The model parameters of the G2++ model, with mean reversion parameters
    [a] and [b], volatility parameters [sigma] and [nu] and correlation
    coefficient [rho]. *)

exception Sdomain
(** This exception is raised by [swaption] if the numerical integration used in the
    price calculation does not converge (indicating that the model parameters are
    malformed). *)

val pricer_of_swaption: today: date -> zc:(date -> float) -> swaption -> (parameters -> float)
val black_price: today: date -> zc: (date -> float) -> swaption -> (float -> float)

val calibrate:
  ?feedback:(string -> unit) ->
    max_global:int ->
      today: date ->
        swaption_quotes: (swaption * float) array ->
          variables:Optimization.optimization_variable array ->
            zc: (date -> float) ->
              unit -> parameters Optimization.calibration_result
