(* ::Package:: *)

(*
  BibergSmoothPipeSolver_v5.wl

  Smooth-wall solver for steady, fully developed,
  stratified turbulent pipe flow with prescribed liquid holdup.

  This version is tailored to co-current positive flow only:
    Ug > 0, Ul > 0, tauI > 0

  Main changes relative to earlier versions:
    - solves in logarithmic variables to enforce positivity
    - uses FindMinimum on the squared scaled residual norm as the primary solver
    - extracts FindMinimum output correctly as replacement rules
    - FindRoot refinement is optional and OFF by default
    - avoids returning the initial guess silently when a good presolve exists
*)

ClearAll["Global`*"];

(* ------------------------------ *)
(* Constants                      *)
(* ------------------------------ *)

kappa = 0.4;
Bconst = 5.5;
small = 10^-10;
epsY = 10^-8;
epsK = 10^-8;

(* ------------------------------ *)
(* Default example case           *)
(* ------------------------------ *)

DefaultBibergCase[] := <|
  "rhoL"   -> 998.2,
  "muL"    -> 1.002*10^-3,
  "rhoG"   -> 1.20,
  "muG"    -> 1.81*10^-5,
  "D"      -> 0.10,
  "alphaL" -> 0.1955,
  "g"      -> 9.81,
  "theta"  -> -2.86 Degree,
  (*"dpdx"   -> -150.0,*)
  "dpdx"   -> -0.0,
  "Ug0"    -> 5.0,
  "Ul0"    -> 1.50,
  "tauI0"  -> 0.5
|>;

(* ------------------------------ *)
(* Geometry from prescribed holdup*)
(* ------------------------------ *)

DeltaFromHoldup[alphaL_?NumericQ] := Module[{delta},
  delta /. FindRoot[
    (delta - Sin[2 delta]/2)/Pi == alphaL,
    {delta, Max[10^-6, Min[Pi - 10^-6, Pi alphaL]]}
  ]
];

PipeGeometry[D_?NumericQ, alphaL_?NumericQ] := Module[
  {deltaL, area, Al, Ag, Si, Swl, Swg, hl, hg},

  deltaL = DeltaFromHoldup[alphaL];
  area = Pi D^2/4;
  Al = alphaL area;
  Ag = (1 - alphaL) area;

  Si = D Sin[deltaL];
  Swl = D deltaL;
  Swg = D (Pi - deltaL);

  hl = Al/Si;
  hg = Ag/Si;

  <|
    "deltaL" -> deltaL,
    "Al" -> Al,
    "Ag" -> Ag,
    "Si" -> Si,
    "Swl" -> Swl,
    "Swg" -> Swg,
    "hl" -> hl,
    "hg" -> hg
  |>
];

(* ------------------------------ *)
(* Smooth-wall friction factor    *)
(* ------------------------------ *)

SmoothFrictionFactor[U_?NumericQ, De_?NumericQ, nu_?NumericQ] := Module[
  {Re, fLam, fTurb},
  Re = Max[U De/nu, 1.0];
  fLam = 64.0/Re;
  fTurb = 1.0/(-1.8 Log10[6.9/Re])^2;
  Max[fLam, fTurb]
];

(* ------------------------------ *)
(* K regularization helpers       *)
(* ------------------------------ *)

RegularizeK[K_?NumericQ, R_?NumericQ] := Module[{k = K, a},
  a = Abs[R]^(5/6);
  k = Which[
    !NumericQ[k], epsK,
    k <= epsK, epsK,
    Abs[k - 1.0] <= epsK, 1.0 - epsK,
    True, k
  ];
  If[Abs[k - a] <= epsK, k = a + epsK];
  k
];

(* ------------------------------ *)
(* Biberg helper functions        *)
(* ------------------------------ *)

DeltaBiberg[Y_?NumericQ, R_?NumericQ, K_?NumericQ] := Module[
  {y = Y, r = R, k = RegularizeK[K, R], a, b, c, d, e, BB, CC, DD, EE},

  Which[
    Abs[r] < 10^-10,
      Log[1 - y] - Log[y + k (1 - y)],

    k < 10^-10,
      Log[1 - y]
      + Sign[r] Sqrt[Abs[r]] Log[y]
      - (1/3) (1 + Sign[r] Sqrt[Abs[r]])
          Log[y^3 + (1 - y)^3 Abs[r]^(5/2)],

    True,
      a  = Abs[r]^(5/6);
      BB = (k^3 + r^3)/(Abs[r]^(5/2) - k^3);
      CC = (r + Sqrt[Abs[r]]) Abs[r]^(1/3)/(3 (k - a));
      DD = -(r + Sqrt[Abs[r]]) (k + 2 a) Abs[r]^(1/3)/
            (6 (k^2 + k a + a^2));
      EE = k (r + Sqrt[Abs[r]]) Abs[r]^(1/3)/
            (Sqrt[3] (k^2 + k a + a^2));

      b = k/(k - 1);
      c = a/(a - 1);
      d = (a + 2 a^2)/(2 (1 + a + a^2));
      e = Sqrt[3] a/(2 (1 + a + a^2));

      Log[1 - y]
      + BB Log[Abs[(y - b)/(1 - b)]]
      + CC Log[Abs[(y - c)/(1 - c)]]
      + DD Log[((y - d)^2 + e^2) (1 + a + a^2)]
      + EE ArcTan[(y - d)/e]
  ]
];

PsiBiberg[R_?NumericQ, K_?NumericQ] := Module[{a, k},
  k = RegularizeK[K, R];
  If[Abs[R] < 10^-10 || k < 10^-10,
    0.0,
    a = Abs[R]^(5/6);
    -k (R + Sqrt[Abs[R]]) Abs[R]^(1/3)/
      (Sqrt[3] (k^2 + k a + a^2)) * ArcTan[(1 + 2 a)/Sqrt[3]]
  ]
];

LambdaFree[KF_: 0.56] := Module[{k = Max[epsK, Min[1 - epsK, KF]]},
  k Log[k]/(1 - k)
];

LambdaBiberg[R_?NumericQ, K_?NumericQ] := Module[{k},
  k = RegularizeK[K, R];
  Which[
    Abs[R] < 10^-10,
      LambdaFree[k],

    True,
      NIntegrate[
        Evaluate[DeltaBiberg[y, R, k]],
        {y, epsY, 1 - epsY},
        Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
        AccuracyGoal -> 8,
        PrecisionGoal -> 8,
        MaxRecursion -> 20
      ]
  ]
];

LambdaP = -(2/27) (2 Sqrt[3] Pi + 9);

ShapeF[R_?NumericQ, K_?NumericQ, KF_: 0.56] := Module[{k},
  k = RegularizeK[K, R];
  (LambdaBiberg[R, k] + PsiBiberg[R, k] - LambdaFree[KF])/
  (LambdaP - LambdaFree[KF])
];

DeffBiberg[A_?NumericQ, Sw_?NumericQ, Si_?NumericQ, R_?NumericQ, K_?NumericQ, KF_: 0.56] := Module[{k},
  k = RegularizeK[K, R];
  4 A/Sw * (Sw/(Sw + Si))^ShapeF[R, k, KF]
];

(* ------------------------------ *)
(* Interface turbulence model     *)
(* ------------------------------ *)

KMin[nu_?NumericQ, rho_?NumericQ, tauI_?NumericQ, h_?NumericQ] :=
  nu Exp[-kappa Bconst]/(Sqrt[Max[tauI, small]/rho] h);

KModel[Ug_?NumericQ, Ul_?NumericQ, geo_Association, case_Association] := Module[
  {nuG, hg, dU},
  nuG = case["muG"]/case["rhoG"];
  hg = geo["hg"];
  dU = Max[Abs[Ug - Ul], 10^-8];
  {
    8 nuG/(dU hg),
    1.0 - epsK
  }
];

(* ------------------------------ *)
(* Residual equations             *)
(* ------------------------------ *)

ResidualsPhysical[
  {Ug_?NumericQ, Ul_?NumericQ, tauI_?NumericQ},
  geo_Association,
  case_Association
] := Module[
  {
    rhoG, rhoL, nuG, nuL, dpdx, grav, theta,
    tauWg, tauWl, Rg, Rl,
    KgRaw, KlRaw, Kg, Kl,
    DeG, DeL, fg, fl, noSlip
  },

  rhoG = case["rhoG"]; rhoL = case["rhoL"];
  nuG = case["muG"]/rhoG;
  nuL = case["muL"]/rhoL;
  dpdx = case["dpdx"];
  grav = case["g"];
  theta = case["theta"];

  tauWg = -(geo["Ag"]/geo["Swg"]) (dpdx + rhoG grav Sin[theta]) - tauI geo["Si"]/geo["Swg"];
  tauWl = -(geo["Al"]/geo["Swl"]) (dpdx + rhoL grav Sin[theta]) + tauI geo["Si"]/geo["Swl"];

  Rg = tauI/tauWg;
  Rl = -tauI/tauWl;

  {KgRaw, KlRaw} = KModel[Ug, Ul, geo, case];

  Kg = RegularizeK[Max[KgRaw, KMin[nuG, rhoG, tauI, geo["hg"]]], Rg];
  Kl = RegularizeK[Max[KlRaw, KMin[nuL, rhoL, tauI, geo["hl"]]], Rl];

  DeG = DeffBiberg[geo["Ag"], geo["Swg"], geo["Si"], Rg, Kg];
  DeL = DeffBiberg[geo["Al"], geo["Swl"], geo["Si"], Rl, Kl];

  fg = SmoothFrictionFactor[Ug, DeG, nuG];
  fl = SmoothFrictionFactor[Ul, DeL, nuL];

  noSlip =
    (Ug - Ul)
    - (
        Sqrt[Abs[tauWg]/rhoG]/kappa
          (LambdaBiberg[Rg, Kg] - DeltaBiberg[epsY, Rg, Kg])
        -
        Sqrt[Abs[tauWl]/rhoL]/kappa
          (LambdaBiberg[Rl, Kl] - DeltaBiberg[epsY, Rl, Kl])
      );

  {
    tauWg - fg rhoG Ug^2/8.0,
    tauWl - fl rhoL Ul^2/8.0,
    noSlip
  }
];

ScaledResidualsPhysical[
  {Ug_?NumericQ, Ul_?NumericQ, tauI_?NumericQ},
  geo_Association,
  case_Association
] := Module[
  {r, rhoG, rhoL, dpdx, grav, theta, tauWgScale, tauWlScale, velScale},

  rhoG = case["rhoG"]; rhoL = case["rhoL"];
  dpdx = case["dpdx"];
  grav = case["g"];
  theta = case["theta"];

  tauWgScale = Max[1.0, Abs[-(geo["Ag"]/geo["Swg"]) (dpdx + rhoG grav Sin[theta])]];
  tauWlScale = Max[1.0, Abs[-(geo["Al"]/geo["Swl"]) (dpdx + rhoL grav Sin[theta])]];
  velScale = Max[1.0, Ug, Ul];

  r = ResidualsPhysical[{Ug, Ul, tauI}, geo, case];
  {r[[1]]/tauWgScale, r[[2]]/tauWlScale, r[[3]]/velScale}
];

ScaledResidualsLog[
  {xg_?NumericQ, xl_?NumericQ, xt_?NumericQ},
  geo_Association,
  case_Association
] := Module[{Ug, Ul, tauI},
  Ug = Exp[xg];
  Ul = Exp[xl];
  tauI = Exp[xt];
  ScaledResidualsPhysical[{Ug, Ul, tauI}, geo, case]
];

ObjectiveLogVars[
  {xg_?NumericQ, xl_?NumericQ, xt_?NumericQ},
  geo_Association,
  case_Association
] := Module[{r},
  r = ScaledResidualsLog[{xg, xl, xt}, geo, case];
  r . r
];

(* ------------------------------ *)
(* Main solver                    *)
(* ------------------------------ *)

Options[SolveBibergSmooth] = {
  WorkingPrecision -> 30,
  AccuracyGoal -> 10,
  PrecisionGoal -> 10,
  MaxIterations -> 100,
  ResidualTolerance -> 10^-6,
  UseFindRootRefinement -> False
};

InitialResiduals[case_Association] := Module[{geo},
  geo = PipeGeometry[case["D"], case["alphaL"]];
  ResidualsPhysical[{case["Ug0"], case["Ul0"], case["tauI0"]}, geo, case]
];

SolveBibergSmooth[case_Association, OptionsPattern[]] := Module[
  {
    geo, xg0, xl0, xt0,
    xg, xl, xt,
    pre, preRules,
    xgSol, xlSol, xtSol,
    UgSol, UlSol, tauISol,
    root, rootRules,
    scaledRes, physRes, solverStatus,
    tauWg, tauWl, Rg, Rl, KgRaw, KlRaw, Kg, Kl, DeG, DeL, fg, fl,
    tol
  },

  geo = PipeGeometry[case["D"], case["alphaL"]];
  tol = OptionValue[ResidualTolerance];

  xg0 = Log[Max[case["Ug0"], 10^-6]];
  xl0 = Log[Max[case["Ul0"], 10^-6]];
  xt0 = Log[Max[case["tauI0"], 10^-6]];

  pre = Quiet @ Check[
    FindMinimum[
      ObjectiveLogVars[{xg, xl, xt}, geo, case],
      {{xg, xg0}, {xl, xl0}, {xt, xt0}},
      WorkingPrecision -> OptionValue[WorkingPrecision],
      AccuracyGoal -> OptionValue[AccuracyGoal],
      PrecisionGoal -> OptionValue[PrecisionGoal],
      MaxIterations -> OptionValue[MaxIterations]
    ],
    $Failed
  ];

  If[pre === $Failed || !ListQ[pre] || Length[pre] < 2 || !ListQ[pre[[2]]],
    Return[<|
      "Ug" -> case["Ug0"],
      "Ul" -> case["Ul0"],
      "tauI" -> case["tauI0"],
      "ResidualsAtSolution" -> InitialResiduals[case],
      "ScaledResidualsAtSolution" -> Missing["NotAvailable"],
      "ResidualNorm" -> Infinity,
      "SolverStatus" -> "FindMinimumFailed",
      "PresolveResult" -> pre,
      "RootObject" -> $Failed
    |>]
  ];

  preRules = pre[[2]];
  xgSol = xg /. preRules;
  xlSol = xl /. preRules;
  xtSol = xt /. preRules;

  UgSol = Exp[xgSol];
  UlSol = Exp[xlSol];
  tauISol = Exp[xtSol];

  scaledRes = ScaledResidualsPhysical[{UgSol, UlSol, tauISol}, geo, case];
  physRes = ResidualsPhysical[{UgSol, UlSol, tauISol}, geo, case];
  solverStatus = "FindMinimumAccepted";

  If[OptionValue[UseFindRootRefinement] && Norm[scaledRes] < 10 tol,
    root = Quiet @ Check[
      FindRoot[
        {
          ScaledResidualsLog[{xg, xl, xt}, geo, case][[1]] == 0,
          ScaledResidualsLog[{xg, xl, xt}, geo, case][[2]] == 0,
          ScaledResidualsLog[{xg, xl, xt}, geo, case][[3]] == 0
        },
        {{xg, xgSol}, {xl, xlSol}, {xt, xtSol}},
        Method -> "Newton",
        WorkingPrecision -> OptionValue[WorkingPrecision],
        AccuracyGoal -> OptionValue[AccuracyGoal],
        PrecisionGoal -> OptionValue[PrecisionGoal],
        MaxIterations -> OptionValue[MaxIterations]
      ],
      $Failed
    ];

    If[root =!= $Failed && ListQ[root],
      rootRules = root;
      xgSol = xg /. rootRules;
      xlSol = xl /. rootRules;
      xtSol = xt /. rootRules;
      UgSol = Exp[xgSol];
      UlSol = Exp[xlSol];
      tauISol = Exp[xtSol];
      scaledRes = ScaledResidualsPhysical[{UgSol, UlSol, tauISol}, geo, case];
      physRes = ResidualsPhysical[{UgSol, UlSol, tauISol}, geo, case];
      solverStatus = "FindMinimumAccepted+FindRootRefined";
      ,
      root = $Failed;
      solverStatus = "FindMinimumAccepted+FindRootFailed";
    ];
    ,
    root = $Failed;
  ];

  tauWg = -(geo["Ag"]/geo["Swg"]) (case["dpdx"] + case["rhoG"] case["g"] Sin[case["theta"]]) - tauISol geo["Si"]/geo["Swg"];
  tauWl = -(geo["Al"]/geo["Swl"]) (case["dpdx"] + case["rhoL"] case["g"] Sin[case["theta"]]) + tauISol geo["Si"]/geo["Swl"];

  Rg = tauISol/tauWg;
  Rl = -tauISol/tauWl;

  {KgRaw, KlRaw} = KModel[UgSol, UlSol, geo, case];
  Kg = RegularizeK[Max[KgRaw, KMin[case["muG"]/case["rhoG"], case["rhoG"], tauISol, geo["hg"]]], Rg];
  Kl = RegularizeK[Max[KlRaw, KMin[case["muL"]/case["rhoL"], case["rhoL"], tauISol, geo["hl"]]], Rl];

  DeG = DeffBiberg[geo["Ag"], geo["Swg"], geo["Si"], Rg, Kg];
  DeL = DeffBiberg[geo["Al"], geo["Swl"], geo["Si"], Rl, Kl];

  fg = SmoothFrictionFactor[UgSol, DeG, case["muG"]/case["rhoG"]];
  fl = SmoothFrictionFactor[UlSol, DeL, case["muL"]/case["rhoL"]];

  <|
    "Ug" -> UgSol,
    "Ul" -> UlSol,
    "tauI" -> tauISol,
    "tauWg" -> tauWg,
    "tauWl" -> tauWl,
    "Rg" -> Rg,
    "Rl" -> Rl,
    "Kg" -> Kg,
    "Kl" -> Kl,
    "DeG" -> DeG,
    "DeL" -> DeL,
    "fg" -> fg,
    "fl" -> fl,
    "Geometry" -> geo,
    "ResidualsAtSolution" -> physRes,
    "ScaledResidualsAtSolution" -> scaledRes,
    "ResidualNorm" -> Norm[scaledRes],
    "SolverStatus" -> solverStatus,
    "PresolveResult" -> pre,
    "RootObject" -> root
  |>
];

RunBibergExample[] := Module[{case, sol},
  case = DefaultBibergCase[];
  sol = SolveBibergSmooth[case];
  <|"Case" -> case, "Solution" -> sol|>
];
