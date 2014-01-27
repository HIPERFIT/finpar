(***************************************************************************)
(*  Copyright (C) 2000-2012 LexiFi SAS. All rights reserved.               *)
(*                                                                         *)
(*  No part of this document may be reproduced or transmitted in any       *)
(*  form or for any purpose without the express permission of LexiFi SAS.  *)
(***************************************************************************)

(* $Id: date.mli 1118 2012-12-05 10:10:06Z afrisch $ *)

type date
val max_date: date
val min_date: date
val of_string: string -> date
val days_between: date -> date -> float
val act_365: date -> date -> float
val add_months: date -> int -> date
val add_years: date -> int -> date
