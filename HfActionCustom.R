# HF-ACTION: multiple imputation for recurrent hospitalizations and mortality.
# Data: rmt::hfaction. Source this file to run the analysis.

HfActionSettings <- list(
  Horizon = 2,
  Intervals = 8,
  Cutpoints = NULL,
  Seed = 29092026,
  CoreFraction = 0.60,
  Chains = 2,
  Warmup = 2000,
  Sampling = 8000,  # Iterations per chain in each assessment or sampling block.
  MaximumWarmupBlocks = 5,
  MaximumSamplingBlocks = 5,
  FixedEffectSd = 5,
  AlphaSd = 5,
  FrailtyShape = 0.5,
  FrailtyRate = 0.05,
  BaselineLogMean = 0,
  BaselineLogSd = 5,
  BlockMoves = TRUE,
  MaximumRhat = 1.05,
  MinimumEss = 400,
  ConfidenceLevel = 0.95,
  ImputationBatchSize = 200,
  OutputFolder = "HfActionCustom",
  InstallMissingPackages = TRUE
)

HfActionProgramPath <- local({
  Paths <- unlist(lapply(sys.frames(), function(Frame) Frame$ofile), use.names = FALSE)
  Arguments <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(Paths)) normalizePath(tail(Paths, 1), winslash = "/", mustWork = TRUE)
  else if (length(Arguments)) normalizePath(sub("^--file=", "", Arguments[1]), winslash = "/", mustWork = TRUE)
  else ""
})


HfActionPlotSettings <- list(
  FontFamily = "sans",
  PdfFallbackFamily = "Helvetica",
  BaseFontSize = 12,
  SamplerFigureWidth = 8.3,
  SamplerFigureHeight = 9.4,
  SamplerLegendSpace = 0.35,
  SamplerYLabelLine = 2.2,
  AxisLabelSize = 12,
  TickLabelSize = 11,
  ScenarioLabelSize = 11,
  LegendSize = 10,
  AnnotationSize = 10,
  PanelLabelSize = 12,
  FigureWidth = 8.3,
  FigureHeight = 3.6,
  MainMargins = c(0.9, 2.25, 0.2, 0.2),
  SupplementSingleWidth = 8.3,
  SupplementSingleHeight = 5.2,
  SupplementWidth = 11.7,
  SupplementHeight = 8.3,
  # Separate model figures: inches; fonts are points at the saved PDF width.
  ModelFigureWidth = 6.7,
  ModelObservedFitHeight = 6.5,
  ModelBaselineHeight = 3.9,
  ModelCoxSnellHeight = 4.2,
  ModelAxisLabelSize = 11,
  ModelTickLabelSize = 10.5,
  ModelLegendSize = 10.5,
  ModelMargins = c(0.8, 0.82, 0.18, 0.12),
  ModelLegendSpace = 0.42,
  ModelXLabelLine = 2.35,
  ModelYLabelLine = 2.8,
  ModelLineWidth = 1.5,
  ModelPointSize = 0.95,
  XLabelLine = 2.6,
  YLabelLine = 3.0,
  LegendBelowInches = 0.64,
  ScenarioLabelOffsetInches = 0.13,
  LineWidth = 1.8,
  AxisLineWidth = 0.8,
  ReferenceLineWidth = 0.9,
  PointSize = 1.1,
  ScatterPointSize = 0.6,
  ScenarioColours = c("#2369A1", "#C47816", "#16836B", "#8B5B9E"),
  ArmColours = c("#2369A1", "#C47816"),
  ReferenceColour = "#777777"
)

# LaTeX settings: font names are standard LaTeX size commands without a backslash.
HfActionTableSettings <- list(
  CountDigits = 3,
  SurvivalDigits = 1,
  TableFont = "small",
  CaptionFont = "small",
  FootnoteFont = "footnotesize",
  ColumnSepPt = 4,
  RowStretch = 1.12
)



CustomSamplerCode <- r"(#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;
// [[Rcpp::plugins(cpp17)]]

struct Evaluation {
  double X, H, D;
};

// Log density relative to the current point. It is a Gaussian log density
// plus a linear term minus positive multiples of exponentials.
struct Conditional {
  double Linear, Precision;
  std::vector<double> Means, Slopes;
  Evaluation Eval(double X) const {
    double H = Linear*X - 0.5*Precision*X*X;
    double D = Linear - Precision*X;
    for (size_t J = 0; J < Means.size(); ++J) {
      if (Means[J] == 0 || Slopes[J] == 0) continue;
      double Y = Slopes[J]*X;
      double V = std::exp(std::log(Means[J]) + Y);
      if (!std::isfinite(V)) return {X, R_NegInf, -std::copysign(R_PosInf, Slopes[J])};
      H -= std::abs(Y) < 0.5 ? Means[J]*std::expm1(Y) : V - Means[J];
      D -= Slopes[J]*V;
    }
    return {X, H, D};
  }
  double Curvature() const {
    double Value = Precision;
    for (size_t J = 0; J < Means.size(); ++J) Value += Means[J]*Slopes[J]*Slopes[J];
    return Value;
  }
};

double OpenUniform() {
  double U;
  do { U = R::runif(0, 1); } while (U <= 0 || U >= 1);
  return U;
}

double LogIntegral(const Evaluation &P, double Left, double Right) {
  if (!std::isfinite(Left)) {
    if (P.D <= 0) stop("The left sampler tail is not integrable.");
    return P.H + P.D*(Right - P.X) - std::log(P.D);
  }
  if (!std::isfinite(Right)) {
    if (P.D >= 0) stop("The right sampler tail is not integrable.");
    return P.H + P.D*(Left - P.X) - std::log(-P.D);
  }
  double Width = Right - Left;
  if (Width <= 0) stop("The sampler envelope has an invalid interval.");
  if (P.D > 0) return P.H + P.D*(Right - P.X) + std::log(-std::expm1(-P.D*Width)) - std::log(P.D);
  if (P.D < 0) return P.H + P.D*(Left - P.X) + std::log(-std::expm1(P.D*Width)) - std::log(-P.D);
  return P.H + std::log(Width);
}

double DrawSegment(double Slope, double Left, double Right) {
  double U = OpenUniform();
  if (!std::isfinite(Left)) return Right + std::log(U)/Slope;
  if (!std::isfinite(Right)) return Left + std::log(U)/Slope;
  double Width = Right - Left;
  if (Slope == 0) return Left + U*Width;
  if (std::abs(Slope*Width) < 0.5) return Left + std::log1p(U*std::expm1(Slope*Width))/Slope;
  if (Slope > 0) return Right + std::log(U + (1 - U)*std::exp(-Slope*Width))/Slope;
  return Left + std::log(U + (1 - U)*std::exp(Slope*Width))/Slope;
}

double DrawConditional(const Conditional &F, double &Evaluations) {
  double Scale = 1/std::sqrt(F.Curvature());
  if (!std::isfinite(Scale) || Scale <= 0) stop("Invalid conditional curvature.");
  Evaluation Centre = F.Eval(0);
  Evaluation Left = F.Eval(-Scale), Right = F.Eval(Scale);
  Evaluations += 3;
  // Gaussian priors guarantee integrable tails. The guard aborts; it never
  // truncates the support or accepts a candidate from an invalid envelope.
  double TailSlope = 1/Scale;
  for (int J = 0; Left.D < TailSlope && J < 100; ++J) {
    Left = F.Eval(Left.X*2); ++Evaluations;
  }
  for (int J = 0; Right.D > -TailSlope && J < 100; ++J) {
    Right = F.Eval(Right.X*2); ++Evaluations;
  }
  if (!std::isfinite(Left.H) || !std::isfinite(Right.H) || Left.D <= 0 || Right.D >= 0) {
    stop("Could not construct the conditional envelope.");
  }
  std::vector<Evaluation> Points = {Left, Centre, Right};
  for (int Attempt = 0; Attempt < 1000; ++Attempt) {
    size_t Size = Points.size();
    std::vector<double> Bounds(Size + 1), LogMass(Size), Mass(Size);
    Bounds[0] = R_NegInf; Bounds[Size] = R_PosInf;
    for (size_t J = 1; J < Size; ++J) {
      const auto &A = Points[J - 1]; const auto &B = Points[J];
      double Difference = A.D - B.D;
      if (Difference <= 0 || !std::isfinite(Difference)) stop("Conditional slopes are not decreasing.");
      double Z = A.X + (B.H - A.H - B.D*(B.X - A.X))/Difference;
      double Tolerance = 1e-9*(1 + std::abs(A.X) + std::abs(B.X));
      if (Z < A.X - Tolerance || Z > B.X + Tolerance || !std::isfinite(Z)) stop("Invalid tangent intersection.");
      Bounds[J] = std::max(A.X, std::min(B.X, Z));
    }
    for (size_t J = 0; J < Size; ++J) LogMass[J] = LogIntegral(Points[J], Bounds[J], Bounds[J + 1]);
    double Max = *std::max_element(LogMass.begin(), LogMass.end()), Total = 0;
    for (size_t J = 0; J < Size; ++J) { Mass[J] = std::exp(LogMass[J] - Max); Total += Mass[J]; }
    if (!std::isfinite(Total) || Total <= 0) stop("Invalid sampler envelope mass.");
    double Target = OpenUniform()*Total, Sum = 0;
    size_t Segment = Size - 1;
    for (size_t J = 0; J < Size; ++J) { Sum += Mass[J]; if (Target <= Sum) { Segment = J; break; } }
    double X = DrawSegment(Points[Segment].D, Bounds[Segment], Bounds[Segment + 1]);
    double Upper = Points[Segment].H + Points[Segment].D*(X - Points[Segment].X);
    Evaluation Candidate = F.Eval(X); ++Evaluations;
    if (std::isnan(Candidate.H) || Candidate.H - Upper > 1e-8*(1 + std::abs(Candidate.H))) {
      stop("The conditional density exceeds its sampler envelope.");
    }
    if (std::log(OpenUniform()) <= Candidate.H - Upper) return X;
    if (!std::isfinite(Candidate.H) || !std::isfinite(Candidate.D)) {
      // An overflowing exponential has zero acceptance probability. Refine
      // the envelope at a finite point between this rejection and its knot.
      for (int Refine = 0; Refine < 100 && (!std::isfinite(Candidate.H) || !std::isfinite(Candidate.D)); ++Refine) {
        X = 0.5*(X + Points[Segment].X);
        Candidate = F.Eval(X); ++Evaluations;
      }
      if (!std::isfinite(Candidate.H) || !std::isfinite(Candidate.D)) stop("Nonfinite conditional evaluation.");
    }
    auto At = std::lower_bound(Points.begin(), Points.end(), X,
      [](const Evaluation &P, double Value) { return P.X < Value; });
    Points.insert(At, Candidate);
  }
  stop("Conditional sampling exceeded its iteration guard.");
  return NA_REAL;
}

// [[Rcpp::export]]
NumericVector SampleConditionalCpp(int Number, double Linear, double Precision,
    NumericVector Means, NumericVector Slopes) {
  if (Number < 1 || Precision <= 0 || Means.size() != Slopes.size()) stop("Invalid conditional input.");
  Conditional F{Linear, Precision, as<std::vector<double>>(Means), as<std::vector<double>>(Slopes)};
  for (double X : F.Means) if (!std::isfinite(X) || X < 0) stop("Invalid conditional weight.");
  NumericVector Out(Number); double Evaluations = 0;
  for (int J = 0; J < Number; ++J) {
    if (J%1000 == 0) checkUserInterrupt();
    Out[J] = DrawConditional(F, Evaluations);
  }
  Out.attr("Evaluations") = Evaluations;
  return Out;
}

// [[Rcpp::export]]
NumericMatrix EvaluateConditionalCpp(NumericVector X, double Linear, double Precision,
    NumericVector Means, NumericVector Slopes) {
  Conditional F{Linear, Precision, as<std::vector<double>>(Means), as<std::vector<double>>(Slopes)};
  NumericMatrix Out(X.size(), 2);
  for (int J = 0; J < X.size(); ++J) { Evaluation E = F.Eval(X[J]); Out(J, 0) = E.H; Out(J, 1) = E.D; }
  return Out;
}

// [[Rcpp::export]]
List GibbsChainCpp(List Data, List Initial, int Warmup, int Sampling, bool BlockMoves = true,
    Nullable<List> Continuation = R_NilValue) {
  int N = as<int>(Data["Participants"]), K = as<int>(Data["Intervals"]);
  int Times = as<int>(Data["Times"]), Cells = as<int>(Data["Cells"]);
  if (N < 2 || K < 1 || Warmup < 0 || Sampling < 1) stop("Invalid sampler dimensions.");
  IntegerVector Id = Data["Id"], Bin = Data["Bin"], CountR = Data["CountR"], CountD = Data["CountD"];
  NumericMatrix X = Data["X"], BaseExposure = Data["BaselineExposure"];
  NumericVector LogExposure = Data["LogExposure"];
  double FixedPrecision = 1/std::pow(as<double>(Data["FixedSd"]), 2);
  double AlphaPrecision = 1/std::pow(as<double>(Data["AlphaSd"]), 2);
  double BasePrecision = 1/std::pow(as<double>(Data["BaselineLogSd"]), 2);
  double BaseMean = as<double>(Data["BaselineLogMean"]);
  double Shape = as<double>(Data["FrailtyShape"]), Rate = as<double>(Data["FrailtyRate"]);
  double Alpha = as<double>(Initial["Alpha"]), Precision = std::exp(as<double>(Initial["LogFrailtyPrecision"]));
  std::vector<double> Coefficients = {as<double>(Initial["BetaR"]), as<double>(Initial["GammaR"]),
    as<double>(Initial["BetaD"]), as<double>(Initial["GammaD"])};
  std::vector<double> R = as<std::vector<double>>(Initial["RecurrenceLogRate"]);
  std::vector<double> D = as<std::vector<double>>(Initial["DeathLogRate"]);
  std::vector<double> B = as<std::vector<double>>(Initial["FrailtyStandard"]);
  std::vector<double> A(N, -1), Z(N, -1), Exposure(N*K, 0), MeanR(N*K, 0), MeanD(N*K, 0);
  std::vector<double> NR(N, 0), ND(N, 0), RCount(K, 0), DCount(K, 0);
  for (double &Value : B) Value /= std::sqrt(Precision);
  for (int J = 0; J < Cells; ++J) {
    int I = Id[J] - 1, H = Bin[J] - 1, P = I*K + H;
    if (I < 0 || I >= N || H < 0 || H >= K || Exposure[P] != 0 ||
        (X(J, 0) != 0 && X(J, 0) != 1) || (X(J, 1) != 0 && X(J, 1) != 1)) stop("Invalid risk cells.");
    if (A[I] >= 0 && (A[I] != X(J, 0) || Z[I] != X(J, 1))) stop("Covariates vary within a participant.");
    Exposure[P] = std::exp(LogExposure[J]); A[I] = X(J, 0); Z[I] = X(J, 1);
    NR[I] += CountR[J]; ND[I] += CountD[J]; RCount[H] += CountR[J]; DCount[H] += CountD[J];
  }
  for (int I = 0; I < N; ++I) if (A[I] < 0 || ND[I] > 1) stop("Invalid participant history.");
  auto RefreshMeans = [&]() {
    for (int I = 0; I < N; ++I) for (int H = 0; H < K; ++H) {
      int P = I*K + H;
      MeanR[P] = Exposure[P]*std::exp(R[H] + Coefficients[0]*A[I] + Coefficients[1]*Z[I] + B[I]);
      MeanD[P] = Exposure[P]*std::exp(D[H] + Coefficients[2]*A[I] + Coefficients[3]*Z[I] + Alpha*B[I]);
      if (!std::isfinite(MeanR[P]) || !std::isfinite(MeanD[P])) stop("Nonfinite event intensity.");
    }
  };
  int Offset = 0;
  if (Continuation.isNotNull()) {
    if (Warmup != 0) stop("Continuation must not repeat initial warmup.");
    List State(Continuation);
    Coefficients = as<std::vector<double>>(State["Coefficients"]);
    R = as<std::vector<double>>(State["RecurrenceLogRate"]);
    D = as<std::vector<double>>(State["DeathLogRate"]);
    B = as<std::vector<double>>(State["Frailty"]);
    Alpha = as<double>(State["Alpha"]);
    Precision = as<double>(State["Precision"]);
    MeanR = as<std::vector<double>>(State["MeanR"]);
    MeanD = as<std::vector<double>>(State["MeanD"]);
    Offset = as<int>(State["Iterations"]);
    if (Coefficients.size() != 4 || R.size() != K || D.size() != K || B.size() != N ||
        MeanR.size() != N*K || MeanD.size() != N*K || Offset < 0 || !std::isfinite(Alpha) ||
        !std::isfinite(Precision) || Precision <= 0) stop("Invalid continuation state.");
  } else {
    RefreshMeans();
  }
  std::vector<double> ContinueR, ContinueD;
  int Variables = 8 + 2*K + 2*Times + N + 7;
  NumericMatrix Draws(Sampling, Variables);
  CharacterVector Names(Variables); int Pos = 0;
  for (const auto &Name : {"BetaR", "GammaR", "BetaD", "GammaD", "Alpha", "LogFrailtyPrecision", "FrailtyVariance", "FrailtyPrecision"}) Names[Pos++] = Name;
  for (const auto &Name : {"RecurrenceLogRate", "DeathLogRate"}) for (int H = 0; H < K; ++H) Names[Pos++] = std::string(Name) + "[" + std::to_string(H + 1) + "]";
  for (const auto &Name : {"RecurrenceCumulative", "DeathCumulative"}) for (int H = 0; H < Times; ++H) Names[Pos++] = std::string(Name) + "[" + std::to_string(H + 1) + "]";
  for (int I = 0; I < N; ++I) Names[Pos++] = "B[" + std::to_string(I + 1) + "]";
  for (int J = 0; J < 5; ++J) Names[Pos++] = "Scores[" + std::to_string(J + 1) + "]";
  Names[Pos++] = "ExpectedTotals[1]"; Names[Pos++] = "ExpectedTotals[2]";
  colnames(Draws) = Names;
  double Evaluations = 0, LargestCacheDifference = 0;
  int Total = Warmup + Sampling;
  for (int Iteration = 0; Iteration < Total; ++Iteration) {
    if (Iteration%100 == 0) checkUserInterrupt();
    for (int Process = 0; Process < 2; ++Process) {
      auto &Means = Process == 0 ? MeanR : MeanD;
      auto &Rates = Process == 0 ? R : D;
      auto &Counts = Process == 0 ? RCount : DCount;
      for (int H = 0; H < K; ++H) {
        double Sum = 0; for (int I = 0; I < N; ++I) Sum += Means[I*K + H];
        Conditional F{Counts[H] - BasePrecision*(Rates[H] - BaseMean), BasePrecision, {Sum}, {1}};
        double Delta = DrawConditional(F, Evaluations); Rates[H] += Delta;
        double Multiplier = std::exp(Delta);
        for (int I = 0; I < N; ++I) Means[I*K + H] *= Multiplier;
      }
      for (int Covariate = 0; Covariate < 2; ++Covariate) {
        auto &Cov = Covariate == 0 ? A : Z;
        auto &CountsById = Process == 0 ? NR : ND;
        double Count = 0, Mean = 0; int J = 2*Process + Covariate;
        for (int I = 0; I < N; ++I) if (Cov[I] == 1) {
          Count += CountsById[I]; for (int H = 0; H < K; ++H) Mean += Means[I*K + H];
        }
        Conditional F{Count - FixedPrecision*Coefficients[J], FixedPrecision, {Mean}, {1}};
        double Delta = DrawConditional(F, Evaluations); Coefficients[J] += Delta;
        double Multiplier = std::exp(Delta);
        for (int I = 0; I < N; ++I) if (Cov[I] == 1) for (int H = 0; H < K; ++H) Means[I*K + H] *= Multiplier;
      }
    }
    Conditional FAlpha{-AlphaPrecision*Alpha, AlphaPrecision, std::vector<double>(N, 0), B};
    for (int I = 0; I < N; ++I) {
      FAlpha.Linear += ND[I]*B[I];
      for (int H = 0; H < K; ++H) FAlpha.Means[I] += MeanD[I*K + H];
    }
    double DeltaAlpha = DrawConditional(FAlpha, Evaluations); Alpha += DeltaAlpha;
    for (int I = 0; I < N; ++I) for (int H = 0; H < K; ++H) MeanD[I*K + H] *= std::exp(DeltaAlpha*B[I]);
    for (int I = 0; I < N; ++I) {
      double SumR = 0, SumD = 0;
      for (int H = 0; H < K; ++H) { SumR += MeanR[I*K + H]; SumD += MeanD[I*K + H]; }
      Conditional F{NR[I] + Alpha*ND[I] - Precision*B[I], Precision, {SumR, SumD}, {1, Alpha}};
      double Delta = DrawConditional(F, Evaluations); B[I] += Delta;
      double MR = std::exp(Delta), MD = std::exp(Alpha*Delta);
      for (int H = 0; H < K; ++H) { MeanR[I*K + H] *= MR; MeanD[I*K + H] *= MD; }
    }
    if (BlockMoves) {
      // Exact Gaussian moves along directions that leave both likelihoods
      // unchanged. These are conditional updates in a linear coordinate system.
      for (int Direction = 0; Direction < 3; ++Direction) {
        double P = 0, L = 0;
        for (int I = 0; I < N; ++I) {
          double Cov = Direction == 0 ? 1 : (Direction == 1 ? A[I] : Z[I]);
          P += Precision*Cov*Cov; L -= Precision*Cov*B[I];
        }
        if (Direction == 0) {
          P += K*BasePrecision*(1 + Alpha*Alpha);
          for (int H = 0; H < K; ++H) L += BasePrecision*((R[H] - BaseMean) + Alpha*(D[H] - BaseMean));
        } else {
          P += FixedPrecision*(1 + Alpha*Alpha);
          L += FixedPrecision*(Coefficients[Direction - 1] + Alpha*Coefficients[Direction + 1]);
        }
        double Shift = R::rnorm(L/P, 1/std::sqrt(P));
        for (int I = 0; I < N; ++I) B[I] += Shift*(Direction == 0 ? 1 : (Direction == 1 ? A[I] : Z[I]));
        if (Direction == 0) for (int H = 0; H < K; ++H) { R[H] -= Shift; D[H] -= Alpha*Shift; }
        else { Coefficients[Direction - 1] -= Shift; Coefficients[Direction + 1] -= Alpha*Shift; }
      }
    }
    double SumSquares = 0; for (double Value : B) SumSquares += Value*Value;
    Precision = R::rgamma(Shape + 0.5*N, 1/(Rate + 0.5*SumSquares));
    if (!std::isfinite(Precision) || Precision <= 0) stop("Invalid sampled frailty precision.");
    if (Iteration == Total - 1) {
      ContinueR = MeanR;
      ContinueD = MeanD;
    }
    if ((Offset + Iteration)%100 == 0 || Iteration == Total - 1) {
      std::vector<double> OldR = MeanR, OldD = MeanD;
      RefreshMeans();
      for (int P = 0; P < N*K; ++P) {
        LargestCacheDifference = std::max(LargestCacheDifference, std::abs(OldR[P] - MeanR[P])/(1 + MeanR[P]));
        LargestCacheDifference = std::max(LargestCacheDifference, std::abs(OldD[P] - MeanD[P])/(1 + MeanD[P]));
      }
      if (LargestCacheDifference > 1e-8) stop("Cached intensities disagree with direct reconstruction.");
    }
    if (Iteration == Total - 1 && (Offset + Iteration)%100 == 0) {
      ContinueR = MeanR;
      ContinueD = MeanD;
    }
    if (Iteration >= Warmup) {
      int Row = Iteration - Warmup; Pos = 0;
      for (double Value : Coefficients) Draws(Row, Pos++) = Value;
      Draws(Row, Pos++) = Alpha; Draws(Row, Pos++) = std::log(Precision);
      Draws(Row, Pos++) = 1/Precision; Draws(Row, Pos++) = Precision;
      for (double Value : R) Draws(Row, Pos++) = Value;
      for (double Value : D) Draws(Row, Pos++) = Value;
      for (int Process = 0; Process < 2; ++Process) for (int T = 0; T < Times; ++T) {
        double Sum = 0; for (int H = 0; H < K; ++H) Sum += BaseExposure(T, H)*std::exp(Process == 0 ? R[H] : D[H]);
        Draws(Row, Pos++) = Sum;
      }
      for (double Value : B) Draws(Row, Pos++) = Value;
      double Scores[5] = {0, 0, 0, 0, 0}, Totals[2] = {0, 0};
      for (int I = 0; I < N; ++I) {
        double SR = 0, SD = 0;
        for (int H = 0; H < K; ++H) { SR += MeanR[I*K + H]; SD += MeanD[I*K + H]; }
        Scores[0] += A[I]*(NR[I] - SR); Scores[1] += Z[I]*(NR[I] - SR);
        Scores[2] += A[I]*(ND[I] - SD); Scores[3] += Z[I]*(ND[I] - SD);
        Scores[4] += B[I]*(ND[I] - SD);
        Totals[0] += SR; Totals[1] += SD;
      }
      for (int J = 0; J < 4; ++J) Draws(Row, Pos++) = Scores[J] - FixedPrecision*Coefficients[J];
      Draws(Row, Pos++) = Scores[4] - AlphaPrecision*Alpha;
      for (double Value : Totals) Draws(Row, Pos++) = Value;
    }
  }
  return List::create(_["Draws"] = Draws, _["ConditionalEvaluations"] = Evaluations,
    _["MaximumCacheDifference"] = LargestCacheDifference,
    _["State"] = List::create(_["Coefficients"] = Coefficients,
      _["RecurrenceLogRate"] = R, _["DeathLogRate"] = D, _["Frailty"] = B,
      _["Alpha"] = Alpha, _["Precision"] = Precision,
      _["MeanR"] = ContinueR, _["MeanD"] = ContinueD, _["Iterations"] = Offset + Total));
}
)"

CompileCustomSampler <- function(Directory, Rebuild = FALSE) {
  dir.create(Directory, recursive = TRUE, showWarnings = FALSE)
  # sourceCpp can update its cache even when loading an existing library.
  # Serialize cache access while the parallel workers start.
  Lock <- filelock::lock(file.path(Directory, "Compilation.lock"), timeout = Inf)
  if (is.null(Lock)) stop("Cannot lock the sampler compilation cache.")
  on.exit(filelock::unlock(Lock), add = TRUE)
  Path <- file.path(Directory, "CustomSampler.cpp")
  Text <- paste0(sub("[\r\n]+$", "", CustomSamplerCode), "\n")
  Existing <- if (file.exists(Path)) readChar(Path, file.info(Path)$size, useBytes = TRUE) else ""
  if (!identical(Existing, Text)) {
    Connection <- file(Path, "wb")
    tryCatch(writeChar(Text, Connection, eos = NULL, useBytes = TRUE), finally = close(Connection))
  }
  Rcpp::sourceCpp(Path, cacheDir = file.path(Directory, "Cache"), rebuild = Rebuild,
    showOutput = Rebuild, verbose = FALSE, env = .GlobalEnv)
  invisible(Path)
}

AuditCustomKernels <- function(D, Initial) {
  C <- CustomComponents(Initial, D)
  Precision <- exp(Initial$LogFrailtyPrecision)
  Rows <- list()
  Add <- function(Name, Linear, PriorPrecision, Means, Slopes, Change) {
    Shifts <- c(-0.4, 0.2, 0.7)
    Actual <- EvaluateConditionalCpp(Shifts, Linear, PriorPrecision, Means, Slopes)[, 1]
    Expected <- vapply(Shifts, function(Shift) CustomDensityR(Change(Shift), D) - CustomDensityR(Initial, D), numeric(1))
    Rows[[length(Rows) + 1]] <<- data.frame(Check = Name, MaximumDifference = max(abs(Actual - Expected)))
  }
  for (J in 1:4) {
    Name <- c("BetaR", "GammaR", "BetaD", "GammaD")[J]
    Covariate <- if (J %in% c(1, 3)) 1 else 2
    Keep <- D$X[, Covariate] == 1
    Counts <- if (J <= 2) D$CountR else D$CountD
    Means <- if (J <= 2) exp(C$EtaR) else exp(C$EtaD)
    Add(Name, sum(Counts[Keep]) - Initial[[Name]]/D$FixedSd^2, 1/D$FixedSd^2,
      sum(Means[Keep]), 1, function(Shift) { P <- Initial; P[[Name]] <- P[[Name]] + Shift; P })
  }
  for (Process in c("RecurrenceLogRate", "DeathLogRate")) {
    Counts <- if (Process == "RecurrenceLogRate") D$CountR else D$CountD
    Means <- if (Process == "RecurrenceLogRate") exp(C$EtaR) else exp(C$EtaD)
    for (K in seq_len(D$Intervals)) {
      Keep <- D$Bin == K
      Add(paste0(Process, K), sum(Counts[Keep]) - (Initial[[Process]][K] - D$BaselineLogMean)/D$BaselineLogSd^2,
        1/D$BaselineLogSd^2, sum(Means[Keep]), 1,
        function(Shift) { P <- Initial; P[[Process]][K] <- P[[Process]][K] + Shift; P })
    }
  }
  Add("Alpha", sum(D$CountD*C$B[D$Id]) - Initial$Alpha/D$AlphaSd^2, 1/D$AlphaSd^2,
    exp(C$EtaD), C$B[D$Id], function(Shift) { P <- Initial; P$Alpha <- P$Alpha + Shift; P })
  for (I in unique(c(1, D$Participants))) {
    Keep <- D$Id == I
    Add(paste0("Frailty", I), sum(D$CountR[Keep]) + Initial$Alpha*sum(D$CountD[Keep]) - Precision*C$B[I],
      Precision, c(sum(exp(C$EtaR[Keep])), sum(exp(C$EtaD[Keep]))), c(1, Initial$Alpha),
      function(Shift) { P <- Initial; P$FrailtyStandard[I] <- P$FrailtyStandard[I] + Shift*sqrt(Precision); P })
  }
  GammaShape <- D$FrailtyShape + D$Participants/2
  GammaRate <- D$FrailtyRate + sum(C$B^2)/2
  Differences <- vapply(c(-0.4, 0.2, 0.7), function(Shift) {
    P <- Initial
    P$LogFrailtyPrecision <- P$LogFrailtyPrecision + Shift
    P$FrailtyStandard <- C$B*exp(P$LogFrailtyPrecision/2)
    Actual <- CustomDensityR(P, D) - CustomDensityR(Initial, D) + (D$Participants/2 - 1)*Shift
    Expected <- dgamma(exp(P$LogFrailtyPrecision), GammaShape, rate = GammaRate, log = TRUE) -
      dgamma(Precision, GammaShape, rate = GammaRate, log = TRUE)
    Actual - Expected
  }, numeric(1))
  Rows[[length(Rows) + 1]] <- data.frame(Check = "Gamma precision conditional", MaximumDifference = max(abs(Differences)))
  PatientRows <- match(seq_len(D$Participants), D$Id)
  for (Direction in 0:2) {
    Cov <- if (Direction == 0) rep(1, D$Participants) else D$X[PatientRows, Direction]
    PP <- Precision*sum(Cov^2)
    LL <- -Precision*sum(Cov*C$B)
    if (Direction == 0) {
      PP <- PP + D$Intervals*(1 + Initial$Alpha^2)/D$BaselineLogSd^2
      LL <- LL + sum(C$R - D$BaselineLogMean + Initial$Alpha*(C$D - D$BaselineLogMean))/D$BaselineLogSd^2
    } else {
      PP <- PP + (1 + Initial$Alpha^2)/D$FixedSd^2
      LL <- LL + (C$Coefficients[Direction] + Initial$Alpha*C$Coefficients[Direction + 2])/D$FixedSd^2
    }
    Differences <- vapply(c(-0.4, 0.2, 0.7), function(Shift) {
      P <- Initial
      P$FrailtyStandard <- (C$B + Shift*Cov)*sqrt(Precision)
      if (Direction == 0) {
        P$RecurrenceLogRate <- C$R - Shift
        P$DeathLogRate <- C$D - Initial$Alpha*Shift
      } else {
        NR <- c("BetaR", "GammaR")[Direction]; ND <- c("BetaD", "GammaD")[Direction]
        P[[NR]] <- P[[NR]] - Shift; P[[ND]] <- P[[ND]] - Initial$Alpha*Shift
      }
      CustomDensityR(P, D) - CustomDensityR(Initial, D) - (LL*Shift - PP*Shift^2/2)
    }, numeric(1))
    Rows[[length(Rows) + 1]] <- data.frame(Check = paste("Gaussian direction", Direction), MaximumDifference = max(abs(Differences)))
  }
  Out <- do.call(rbind, Rows)
  if (any(!is.finite(Out$MaximumDifference)) || any(Out$MaximumDifference > 1e-7)) stop("A custom conditional differs from the joint density.")
  Out
}



SetSimulationSeed <- function(Seed) {
  set.seed(Seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}


SimulationHash <- function(Object) {
  Path <- tempfile("SimulationHash")
  on.exit(unlink(Path), add = TRUE)
  saveRDS(Object, Path, compress = FALSE, version = 2)
  unname(tools::md5sum(Path))
}


SaveSimulationRds <- function(Object, Path) {
  Partial <- paste0(Path, ".part")
  saveRDS(Object, Partial, compress = "gzip", version = 2)
  Backup <- paste0(Path, ".previous")
  if (file.exists(Backup)) unlink(Backup)
  if (file.exists(Path) && !file.rename(Path, Backup)) stop("Cannot replace a checkpoint. Close files using it.")
  if (!file.rename(Partial, Path)) {
    if (file.exists(Backup)) file.rename(Backup, Path)
    stop("Cannot save the checkpoint.")
  }
  if (file.exists(Backup)) unlink(Backup)
  invisible(Path)
}


BuildCustomData <- function(Data, Scenario, Settings = HfActionSettings) {
  K <- length(Scenario$Cutpoints) - 1
  Exposure <- t(vapply(Scenario$BaselineTimes, function(Time) {
    pmax(0, pmin(Time, tail(Scenario$Cutpoints, -1)) - head(Scenario$Cutpoints, -1))
  }, numeric(K)))
  list(Participants = nrow(Data$Patients), Cells = nrow(Data$Cells), Intervals = K,
    Times = length(Scenario$BaselineTimes), Id = as.integer(Data$Cells$Id),
    Bin = as.integer(Data$Cells$Interval), CountR = as.integer(Data$Cells$Recurrence),
    CountD = as.integer(Data$Cells$Death), X = as.matrix(Data$Cells[, c("Treatment", "Age60")]),
    LogExposure = log(Data$Cells$Exposure), BaselineExposure = Exposure,
    BaselineLogMean = Settings$BaselineLogMean, BaselineLogSd = Settings$BaselineLogSd,
    FixedSd = Settings$FixedEffectSd, AlphaSd = Settings$AlphaSd,
    FrailtyShape = Settings$FrailtyShape, FrailtyRate = Settings$FrailtyRate)
}


CustomComponents <- function(P, D) {
  B <- P$FrailtyStandard*exp(-P$LogFrailtyPrecision/2)
  Coefficients <- c(P$BetaR, P$GammaR, P$BetaD, P$GammaD)
  EtaR <- D$LogExposure + P$RecurrenceLogRate[D$Bin] + B[D$Id] +
    as.vector(D$X %*% Coefficients[1:2])
  EtaD <- D$LogExposure + P$DeathLogRate[D$Bin] + P$Alpha*B[D$Id] +
    as.vector(D$X %*% Coefficients[3:4])
  list(R = P$RecurrenceLogRate, D = P$DeathLogRate, B = B,
    EtaR = EtaR, EtaD = EtaD, Coefficients = Coefficients)
}


CustomDensityR <- function(P, D) {
  C <- CustomComponents(P, D)
  LogPrecision <- P$LogFrailtyPrecision
  Value <- sum(stats::dpois(D$CountR, exp(C$EtaR), log = TRUE)) +
    sum(stats::dpois(D$CountD, exp(C$EtaD), log = TRUE)) +
    sum(stats::dnorm(C$Coefficients, 0, D$FixedSd, log = TRUE)) +
    stats::dnorm(P$Alpha, 0, D$AlphaSd, log = TRUE) +
    sum(stats::dnorm(c(C$R, C$D), D$BaselineLogMean, D$BaselineLogSd, log = TRUE))
  # Evaluate the density in precision and actual frailties, then add the
  # Jacobians for log precision and the standardized Gaussian frailties.
  Value + stats::dgamma(exp(LogPrecision), D$FrailtyShape, D$FrailtyRate, log = TRUE) +
    LogPrecision + sum(stats::dnorm(C$B, 0, exp(-LogPrecision/2), log = TRUE)) -
    D$Participants*LogPrecision/2
}


GetCustomInitial <- function(Data, D) {
  RateR <- max(sum(Data$Cells$Recurrence)/sum(Data$Patients$Stop), 0.01)
  RateD <- max(sum(Data$Patients$Death)/sum(Data$Patients$Stop), 0.001)
  list(BetaR = 0, GammaR = 0, BetaD = 0, GammaD = 0, Alpha = 0.3,
    LogFrailtyPrecision = 0, RecurrenceLogRate = rep(log(RateR), D$Intervals),
    DeathLogRate = rep(log(RateD), D$Intervals), FrailtyStandard = rep(0, D$Participants))
}


SummarizeCustomMatrix <- function(X) {
  Q <- stats::quantile(as.numeric(X), c(0.025, 0.5, 0.975), names = FALSE)
  Qse <- posterior::mcse_quantile(X, probs = c(0.025, 0.975), names = FALSE)
  c(Mean = mean(X), Sd = stats::sd(as.numeric(X)), Lower = Q[1], Median = Q[2], Upper = Q[3],
    MeanMonteCarloSe = posterior::mcse_mean(X), SdMonteCarloSe = posterior::mcse_sd(X),
    LowerMonteCarloSe = Qse[1], UpperMonteCarloSe = Qse[2],
    Rhat = posterior::rhat(X), BulkEss = posterior::ess_bulk(X), TailEss = posterior::ess_tail(X))
}


BindSimulationRows <- function(Rows) {
  Rows <- Filter(function(X) is.data.frame(X) && nrow(X) > 0, Rows)
  if (!length(Rows)) return(data.frame())
  Columns <- unique(unlist(lapply(Rows, names), use.names = FALSE))
  Rows <- lapply(Rows, function(X) {
    for (Name in setdiff(Columns, names(X))) X[[Name]] <- NA
    X[Columns]
  })
  Output <- do.call(rbind, Rows)
  rownames(Output) <- NULL
  Output
}


EnsureCustomPackages <- function(Settings) {
  Packages <- c("Rcpp", "posterior", "filelock")
  Missing <- Packages[!vapply(Packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(Missing) && Settings$InstallMissingPackages) {
    message("Installing: ", paste(Missing, collapse = ", "))
    utils::install.packages(Missing, repos = "https://cloud.r-project.org")
  }
  Missing <- Packages[!vapply(Packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(Missing)) stop(paste("Install these packages:", paste(Missing, collapse = ", ")), call. = FALSE)
  list(Versions = c(vapply(Packages, function(P) as.character(utils::packageVersion(P)), character(1)),
    R = as.character(getRversion()), Platform = R.version$platform))
}


JitterCustomInitial <- function(P, Seed, Chain) {
  SetSimulationSeed(Seed)
  for (Name in c("BetaR", "GammaR", "BetaD", "GammaD")) P[[Name]] <- stats::rnorm(1, 0, 0.35)
  P$Alpha <- c(-0.5, 1.5)[1 + (Chain - 1)%%2] + stats::rnorm(1, 0, 0.2)
  P$LogFrailtyPrecision <- stats::runif(1, -0.7, 0.7)
  P$RecurrenceLogRate <- P$RecurrenceLogRate + stats::rnorm(length(P$RecurrenceLogRate), 0, 0.35)
  P$DeathLogRate <- P$DeathLogRate + stats::rnorm(length(P$DeathLogRate), 0, 0.35)
  P$FrailtyStandard <- stats::rnorm(length(P$FrailtyStandard), 0, 0.35)
  P
}


RunCustomChain <- function(Job) {
  if (is.null(Job$Continuation)) SetSimulationSeed(Job$Seed) else {
    assign(".Random.seed", Job$Continuation$RngState, envir = .GlobalEnv)
  }
  Start <- proc.time()[3]
  Warnings <- character()
  Answer <- tryCatch(withCallingHandlers({
    Result <- GibbsChainCpp(Job$Data, Job$Initial, Job$Warmup, Job$Sampling, Job$BlockMoves,
      if (is.null(Job$Continuation)) NULL else Job$Continuation$State)
    Result$RngState <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    Result$Succeeded <- TRUE
    Result
  }, warning = function(W) {
    Warnings <<- c(Warnings, conditionMessage(W))
    invokeRestart("muffleWarning")
  }), error = function(E) list(Succeeded = FALSE, Error = conditionMessage(E)))
  Answer$Chain <- Job$Chain
  Answer$Seed <- Job$Seed
  Answer$Seconds <- unname(proc.time()[3] - Start)
  Answer$Warnings <- unique(Warnings)
  Answer$WorkerPid <- Sys.getpid()
  Answer
}


AuditReturnedStates <- function(Draws, D) {
  Largest <- 0
  for (Chain in seq_len(dim(Draws)[2])) for (Iteration in unique(c(1, dim(Draws)[1]))) {
    X <- Draws[Iteration, Chain, ]
    P <- as.list(X[c("BetaR", "GammaR", "BetaD", "GammaD", "Alpha", "LogFrailtyPrecision")])
    P$RecurrenceLogRate <- X[paste0("RecurrenceLogRate[", seq_len(D$Intervals), "]")]
    P$DeathLogRate <- X[paste0("DeathLogRate[", seq_len(D$Intervals), "]")]
    B <- X[paste0("B[", seq_len(D$Participants), "]")]
    P$FrailtyStandard <- B*exp(P$LogFrailtyPrecision/2)
    C <- CustomComponents(P, D)
    ResidualR <- D$CountR - exp(C$EtaR)
    ResidualD <- D$CountD - exp(C$EtaD)
    Scores <- c(as.vector(crossprod(D$X, ResidualR)) - C$Coefficients[1:2]/D$FixedSd^2,
      as.vector(crossprod(D$X, ResidualD)) - C$Coefficients[3:4]/D$FixedSd^2,
      sum(B[D$Id]*ResidualD) - P$Alpha/D$AlphaSd^2)
    Expected <- c(exp(-P$LogFrailtyPrecision), exp(P$LogFrailtyPrecision),
      as.vector(D$BaselineExposure %*% exp(P$RecurrenceLogRate)),
      as.vector(D$BaselineExposure %*% exp(P$DeathLogRate)), Scores,
      sum(exp(C$EtaR)), sum(exp(C$EtaD)))
    Names <- c("FrailtyVariance", "FrailtyPrecision", paste0("RecurrenceCumulative[", seq_len(D$Times), "]"),
      paste0("DeathCumulative[", seq_len(D$Times), "]"), paste0("Scores[", 1:5, "]"), paste0("ExpectedTotals[", 1:2, "]"))
    Difference <- max(abs(X[Names] - Expected)/(1 + abs(Expected)))
    if (!is.finite(Difference) || Difference > 1e-8) stop("Returned sampler summaries disagree with reconstruction.")
    Largest <- max(Largest, Difference)
  }
  Largest
}


GetCustomDiagnostics <- function(Draws, Settings) {
  if (length(dim(Draws)) != 3 || dim(Draws)[2] != Settings$Chains ||
      any(!is.finite(Draws))) stop("MCMC draws are incomplete or nonfinite.")
  Variables <- do.call(rbind, lapply(dimnames(Draws)[[3]], function(Name) {
    X <- Draws[, , Name]
    data.frame(Variable = Name, Rhat = posterior::rhat(X),
      BulkEss = posterior::ess_bulk(X), TailEss = posterior::ess_tail(X))
  }))
  Variables$Review <- !is.finite(Variables$Rhat) | Variables$Rhat > Settings$MaximumRhat |
    !is.finite(Variables$BulkEss) | Variables$BulkEss < Settings$MinimumEss |
    !is.finite(Variables$TailEss) | Variables$TailEss < Settings$MinimumEss
  Variables
}


GetCustomDiagnosticStatus <- function(Variables, Settings) {
  RhatPassed <- nrow(Variables) > 0 && all(is.finite(Variables$Rhat) & Variables$Rhat <= Settings$MaximumRhat)
  EssPassed <- nrow(Variables) > 0 && all(is.finite(Variables$BulkEss) & Variables$BulkEss >= Settings$MinimumEss &
    is.finite(Variables$TailEss) & Variables$TailEss >= Settings$MinimumEss)
  list(RhatPassed = RhatPassed, EssPassed = EssPassed)
}


GetCustomBlockRecord <- function(Variables, Settings, Block, Stage, StageBlock, DrawsAssessed) {
  Status <- GetCustomDiagnosticStatus(Variables, Settings)
  FiniteMax <- function(X) if (any(is.finite(X))) max(X[is.finite(X)]) else NA_real_
  FiniteMin <- function(X) if (any(is.finite(X))) min(X[is.finite(X)]) else NA_real_
  data.frame(Block = Block, Stage = Stage, StageBlock = StageBlock,
    AddedPerChain = Settings$Sampling, AssessedPerChain = DrawsAssessed,
    MaximumRhat = FiniteMax(Variables$Rhat), MinimumBulkEss = FiniteMin(Variables$BulkEss),
    MinimumTailEss = FiniteMin(Variables$TailEss), RhatPassed = Status$RhatPassed,
    EssPassed = Status$EssPassed, Unavailable = sum(!is.finite(Variables$Rhat) |
      !is.finite(Variables$BulkEss) | !is.finite(Variables$TailEss)), UsedInSummary = FALSE)
}


CollectCustomDraws <- function(Runs, Settings) {
  Failed <- !vapply(Runs, function(X) isTRUE(X$Succeeded), logical(1))
  if (any(Failed)) stop(paste("Chain", Runs[[which(Failed)[1]]]$Chain, "failed:", Runs[[which(Failed)[1]]]$Error))
  if (length(Runs) != Settings$Chains) stop("Some chains are missing.")
  Names <- colnames(Runs[[1]]$Draws)
  Draws <- array(NA_real_, c(Settings$Sampling, Settings$Chains, length(Names)),
    dimnames = list(NULL, NULL, Names))
  for (J in seq_len(Settings$Chains)) {
    if (nrow(Runs[[J]]$Draws) != Settings$Sampling || !identical(colnames(Runs[[J]]$Draws), Names)) {
      stop("A chain returned incomplete draws.")
    }
    Draws[, J, ] <- Runs[[J]]$Draws
  }
  if (any(!is.finite(Draws))) stop("A chain returned nonfinite draws.")
  Draws
}


AppendCustomDraws <- function(Previous, Additional) {
  N <- dim(Previous)[1]
  Added <- dim(Additional)[1]
  if (!identical(dim(Previous)[2:3], dim(Additional)[2:3]) ||
      !identical(dimnames(Previous)[[3]], dimnames(Additional)[[3]])) stop("Continuation variables differ.")
  Out <- array(NA_real_, c(N + Added, dim(Previous)[2:3]),
    dimnames = list(NULL, NULL, dimnames(Previous)[[3]]))
  Out[seq_len(N), , ] <- Previous
  Out[N + seq_len(Added), , ] <- Additional
  Out
}


RunCustomAutoFit <- function(Jobs, Settings) {
  State <- get("CustomWorkerState", envir = .GlobalEnv)
  History <- list()
  Totals <- lapply(Jobs, function(Job) list(Chain = Job$Chain, Seed = Job$Seed,
    Seconds = 0, ConditionalEvaluations = 0, MaximumCacheDifference = 0,
    WorkerPid = NA_real_))
  SamplingSeconds <- 0
  DiagnosticSeconds <- 0
  ExtractionSeconds <- 0
  Block <- 0
  Retained <- NULL
  FinalVariables <- NULL
  WarmupEstablished <- FALSE
  WarmupBlocks <- 0
  SamplingBlocks <- 0
  NextBlock <- function() {
    Start <- proc.time()[3]
    Runs <- if (is.null(State$Cluster)) lapply(Jobs, RunCustomChain) else parallel::parLapply(State$Cluster, Jobs, RunCustomChain)
    SamplingSeconds <<- SamplingSeconds + unname(proc.time()[3] - Start)
    for (Run in Runs) for (Message in Run$Warnings) warning(Message, call. = FALSE)
    Start <- proc.time()[3]
    Draws <- CollectCustomDraws(Runs, Settings)
    for (J in seq_along(Jobs)) {
      Jobs[[J]]$Continuation <<- list(State = Runs[[J]]$State, RngState = Runs[[J]]$RngState)
      Jobs[[J]]$Warmup <<- 0
      Totals[[J]]$Seconds <<- Totals[[J]]$Seconds + Runs[[J]]$Seconds
      Totals[[J]]$ConditionalEvaluations <<- Totals[[J]]$ConditionalEvaluations + Runs[[J]]$ConditionalEvaluations
      Totals[[J]]$MaximumCacheDifference <<- max(Totals[[J]]$MaximumCacheDifference, Runs[[J]]$MaximumCacheDifference)
      Totals[[J]]$WorkerPid <<- Runs[[J]]$WorkerPid
      Totals[[J]]$TotalIterations <<- Runs[[J]]$State$Iterations
    }
    Block <<- Block + 1
    ExtractionSeconds <<- ExtractionSeconds + unname(proc.time()[3] - Start)
    Draws
  }
  Diagnose <- function(Draws) {
    Start <- proc.time()[3]
    Out <- GetCustomDiagnostics(Draws, Settings)
    DiagnosticSeconds <<- DiagnosticSeconds + unname(proc.time()[3] - Start)
    Out
  }
  for (Attempt in seq_len(Settings$MaximumWarmupBlocks)) {
    Latest <- NextBlock()
    Variables <- Diagnose(Latest)
    History[[Block]] <- GetCustomBlockRecord(Variables, Settings, Block,
      "WarmupAssessment", Attempt, dim(Latest)[1])
    WarmupBlocks <- Attempt
    ReportApplicationBlock(History[[Block]])
    Status <- GetCustomDiagnosticStatus(Variables, Settings)
    if (Status$RhatPassed) {
      WarmupEstablished <- TRUE
      Retained <- Latest
      FinalVariables <- Variables
      SamplingBlocks <- 1
      break
    }
    # Provisional blocks are replaced only before a warmup boundary is established.
    Retained <- Latest
    FinalVariables <- Variables
  }
  FirstRetainedBlock <- Block
  if (WarmupEstablished) {
    Status <- GetCustomDiagnosticStatus(FinalVariables, Settings)
    while (!(Status$RhatPassed && Status$EssPassed) && SamplingBlocks < Settings$MaximumSamplingBlocks) {
      Additional <- NextBlock()
      Start <- proc.time()[3]
      Retained <- AppendCustomDraws(Retained, Additional)
      ExtractionSeconds <- ExtractionSeconds + unname(proc.time()[3] - Start)
      SamplingBlocks <- SamplingBlocks + 1
      FinalVariables <- Diagnose(Retained)
      History[[Block]] <- GetCustomBlockRecord(FinalVariables, Settings, Block,
        "RetainedSampling", SamplingBlocks, dim(Retained)[1])
      ReportApplicationBlock(History[[Block]])
      Status <- GetCustomDiagnosticStatus(FinalVariables, Settings)
    }
  } else {
    SamplingBlocks <- 1
  }
  Status <- GetCustomDiagnosticStatus(FinalVariables, Settings)
  StopReason <- if (!WarmupEstablished) "WarmupLimit" else if (Status$RhatPassed && Status$EssPassed) "CriteriaMet" else "SamplingLimit"
  History <- do.call(rbind, History)
  History$UsedInSummary <- History$Block >= FirstRetainedBlock
  SummaryWarmup <- Settings$Warmup + (FirstRetainedBlock - 1)*Settings$Sampling
  for (J in seq_along(Totals)) {
    Totals[[J]]$Retained <- dim(Retained)[1]
    Totals[[J]]$WarmupIterations <- SummaryWarmup
    Totals[[J]]$Blocks <- Block
  }
  list(Draws = Retained, Runs = Totals, VariableDiagnostics = FinalVariables,
    BlockDiagnostics = History, WarmupEstablished = WarmupEstablished,
    WarmupBlocks = WarmupBlocks, SamplingBlocks = SamplingBlocks,
    WarmupIterations = SummaryWarmup, DrawsPerChain = dim(Retained)[1],
    TotalBlocks = Block, StopReason = StopReason, RhatPassed = Status$RhatPassed,
    EssPassed = Status$EssPassed, SamplingSeconds = SamplingSeconds,
    ExtractionSeconds = ExtractionSeconds, DiagnosticSeconds = DiagnosticSeconds)
}


CreateOutputDirectory <- function(OutputDir, Prefix) {
  if (!is.character(OutputDir) || length(OutputDir) != 1 || is.na(OutputDir) || !nzchar(OutputDir)) {
    stop("OutputDir must be one nonempty folder path.", call. = FALSE)
  }
  if (!dir.exists(OutputDir) && !dir.create(OutputDir, recursive = TRUE, showWarnings = FALSE) &&
      !dir.exists(OutputDir)) stop("Could not create the output location.", call. = FALSE)
  Index <- 1
  repeat {
    Directory <- file.path(OutputDir, sprintf("%s%03d", Prefix, Index))
    if (!file.exists(Directory)) {
      if (dir.create(Directory, showWarnings = FALSE)) {
        return(normalizePath(Directory, winslash = "/", mustWork = TRUE))
      }
      # Another process may have claimed this folder; continue without overwriting it.
      if (!file.exists(Directory)) stop("Could not create the output folder.", call. = FALSE)
    }
    Index <- Index + 1
  }
}


GetIntervalExposure <- function(Start, Stop, Cutpoints) {
  if (length(Start) != length(Stop) || any(!is.finite(c(Start, Stop))) ||
      any(Start < 0) || any(Stop < Start) || any(Stop > tail(Cutpoints, 1))) {
    stop("Risk interval times are invalid.", call. = FALSE)
  }
  Exposure <- vapply(seq_len(length(Cutpoints) - 1), function(Index) {
    pmax(0, pmin(Stop, Cutpoints[Index + 1]) - pmax(Start, Cutpoints[Index]))
  }, numeric(length(Start)))
  matrix(Exposure, nrow = length(Start), ncol = length(Cutpoints) - 1)
}


ReadData <- function(Data = NULL) {
  if (is.null(Data)) {
    if (!requireNamespace("rmt", quietly = TRUE)) {
      message("Installing rmt from CRAN to load the public dataset...")
      utils::install.packages("rmt", repos = "https://cloud.r-project.org")
    }
    if (!requireNamespace("rmt", quietly = TRUE)) stop("rmt could not be loaded. Check the installation message above.")
    Holder <- new.env(parent = emptyenv())
    utils::data("hfaction", package = "rmt", envir = Holder)
    Data <- Holder$hfaction
  }
  Input <- Data
  Required <- c("patid", "time", "status", "trt_ab", "age60")
  if (!is.data.frame(Input) || !all(Required %in% names(Input))) {
    stop("The dataset must contain the original rmt fields.")
  }
  Raw <- Input
  D <- data.frame(OriginalId = as.character(Input$patid), Time = Input$time,
                  Status = Input$status, Treatment = Input$trt_ab, Age60 = Input$age60,
                  stringsAsFactors = FALSE)
  if (anyNA(D) || any(!nzchar(D$OriginalId)) || anyDuplicated(D) ||
      !all(vapply(D[-1], is.numeric, logical(1))) || any(!is.finite(as.matrix(D[-1]))) ||
      any(D$Time < 0) || any(!D$Status %in% 0:2) ||
      any(!D$Treatment %in% 0:1) || any(!D$Age60 %in% 0:1)) {
    stop("Missing, duplicated or invalid input records; review the data before fitting.")
  }
  Names <- sort(unique(D$OriginalId))
  D$Id <- match(D$OriginalId, Names)
  D <- D[order(D$Id, D$Time, D$Status == 0), ]
  P <- do.call(rbind, lapply(split(D, D$Id), function(X) {
    End <- X[X$Status %in% c(0, 2), ]
    if (nrow(End) != 1 || End$Time <= 0 || any(X$Time > End$Time) ||
        length(unique(X$Treatment)) != 1 || length(unique(X$Age60)) != 1) {
      stop(paste("Invalid endpoint or changing baseline covariates for", X$OriginalId[1]))
    }
    data.frame(Id = X$Id[1], OriginalId = X$OriginalId[1], Treatment = X$Treatment[1],
               Age60 = X$Age60[1], Start = 0, Stop = End$Time, Death = as.numeric(End$Status == 2),
               Count = sum(X$Status == 1), stringsAsFactors = FALSE)
  }))
  P <- P[order(P$Id), ]
  rownames(P) <- NULL
  if (length(unique(P$Treatment)) != 2 || length(unique(P$Age60)) != 2 ||
      sum(P$Death) < 2 || sum(P$Count) < 2) stop("Insufficient covariate or event variation.")
  Boundary <- D[D$Status == 1 & (D$Time == 0 | D$Time == P$Stop[match(D$Id, P$Id)]), ]
  if (nrow(Boundary)) Boundary$Convention <- ifelse(Boundary$Time == 0,
                                                    "Retained in first baseline interval; raw time unchanged", "Included at observation endpoint")
  list(RawData = Raw, EventData = D, Patients = P, Boundary = Boundary)
}


MakeRiskCells <- function(Prepared, Cutpoints) {
  P <- Prepared$Patients
  D <- Prepared$EventData
  Exposure <- GetIntervalExposure(P$Start, P$Stop, Cutpoints)
  Index <- which(Exposure > 0, arr.ind = TRUE)
  Cells <- data.frame(Id = P$Id[Index[, 1]], Interval = Index[, 2],
                      Exposure = Exposure[Index], Recurrence = 0, Death = 0)
  Cells <- Cells[order(Cells$Id, Cells$Interval), ]
  Keys <- paste(Cells$Id, Cells$Interval, sep = ":")
  for (Endpoint in 1:2) {
    Events <- D[D$Status == Endpoint, ]
    Bins <- as.integer(cut(Events$Time, breaks = Cutpoints, right = TRUE, include.lowest = TRUE))
    Where <- match(paste(Events$Id, Bins, sep = ":"), Keys)
    if (anyNA(Where)) stop("An event has no positive-exposure risk cell.")
    Cells[[c("Recurrence", "Death")[Endpoint]]] <- tabulate(Where, nbins = nrow(Cells))
  }
  if (sum(Cells$Recurrence) != sum(P$Count) || sum(Cells$Death) != sum(P$Death) ||
      abs(sum(Cells$Exposure) - sum(P$Stop)) > 1e-9) stop("Counts or follow-up were lost.")
  Cells$Treatment <- P$Treatment[match(Cells$Id, P$Id)]
  Cells$Age60 <- P$Age60[match(Cells$Id, P$Id)]
  Cells
}


GetNelsonAalen <- function(Time, Event) {
  if (length(Time) != length(Event) || any(!is.finite(Time)) || any(Time < 0) ||
      any(!Event %in% 0:1)) stop("Invalid residual times or event indicators.")
  Times <- sort(unique(Time[Event == 1]))
  if (!length(Times)) return(data.frame(Residual = 0, AtRisk = length(Time), Events = 0, CumulativeHazard = 0))
  Risk <- vapply(Times, function(T) sum(Time >= T), numeric(1))
  Deaths <- vapply(Times, function(T) sum(Time == T & Event == 1), numeric(1))
  # Event and censoring ties: everyone with Time >= T is in the risk set.
  data.frame(Residual = c(0, Times), AtRisk = c(length(Time), Risk),
             Events = c(0, Deaths), CumulativeHazard = c(0, cumsum(Deaths/Risk)))
}


GetMiScenarios <- function() {
  data.frame(Scenario = c("MAR", "J2RHospitalization", "J2RDeath", "J2RBoth"),
             RecurrenceBorrowing = c(0, 1, 0, 1), DeathBorrowing = c(0, 0, 1, 1),
             Description = c("Assigned treatment for both outcomes after censoring",
                             "Reference hospitalization intensity; assigned-treatment death hazard",
                             "Assigned-treatment hospitalization intensity; reference death hazard",
                             "Reference treatment components for both outcomes after censoring"), stringsAsFactors = FALSE)
}


CheckMiScenarios <- function(Scenarios) {
  if (!is.data.frame(Scenarios) || nrow(Scenarios) < 1 ||
      !all(c("Scenario", "RecurrenceBorrowing") %in% names(Scenarios))) stop("Invalid sensitivity scenarios.")
  if (!"DeathBorrowing" %in% names(Scenarios)) stop("Specify both borrowing weights.")
  if (anyNA(Scenarios) || anyDuplicated(Scenarios$Scenario) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9]*$", Scenarios$Scenario))) stop("Invalid scenario names.")
  for (Name in c("RecurrenceBorrowing", "DeathBorrowing")) {
    X <- Scenarios[[Name]]
    if (!is.numeric(X) || any(!is.finite(X)) || any(X < 0 | X > 1)) {
      stop("Borrowing weights must be between zero and one.")
    }
  }
  if (sum(Scenarios$Scenario == "MAR" & Scenarios$RecurrenceBorrowing == 0 &
          Scenarios$DeathBorrowing == 0) != 1) stop("Include MAR with both borrowing weights zero.")
  Scenarios
}


PrepareMiTarget <- function(Prepared, Horizon) {
  if (length(Horizon) != 1 || !is.finite(Horizon) || Horizon <= 0 ||
      Horizon > max(Prepared$Patients$Stop)) stop("Horizon must be positive and within observed follow-up.")
  P <- Prepared$Patients
  D <- Prepared$EventData
  Events <- D[D$Status == 1 & D$Time <= Horizon, ]
  P$ObservedCount <- tabulate(match(Events$Id, P$Id), nbins = nrow(P))
  P$ObservedEnd <- pmin(P$Stop, Horizon)
  P$KnownDeath <- as.numeric(P$Death == 1 & P$Stop <= Horizon)
  P$NeedsImputation <- P$Death == 0 & P$Stop < Horizon
  P$Horizon <- Horizon
  P$OutcomeStatus <- ifelse(P$NeedsImputation, "Censored before horizon",
                            ifelse(P$KnownDeath == 1, "Observed death by horizon", "Observed through horizon"))
  P
}


CheckMiDraws <- function(Draws, Target, Cutpoints) {
  M <- nrow(Draws$Fixed)
  N <- nrow(Target)
  K <- length(Cutpoints) - 1
  Required <- c("Treatment_S1", "Age60_S1", "Treatment_S2", "Age60_S2")
  if (length(M) != 1 || M < 2 || !all(Required %in% colnames(Draws$Fixed))) stop("Missing joint fixed-effect draws.")
  for (Name in c("RecurrenceBaseline", "DeathBaseline", "RecurrenceFrailty", "DeathFrailty")) {
    Rows <- if (grepl("Baseline", Name)) K else N
    X <- Draws[[Name]]
    if (!is.matrix(X) || !identical(as.numeric(dim(X)), as.numeric(c(Rows, M))) || any(!is.finite(X))) {
      stop(paste("Invalid joint draws for", Name))
    }
  }
  if (any(!is.finite(Draws$Fixed[, Required])) || anyNA(Target) ||
      any(Target$ObservedCount < 0 | Target$ObservedCount != floor(Target$ObservedCount)) ||
      any(Target$ObservedEnd <= 0 | Target$ObservedEnd > Target$Horizon) ||
      any(Target$KnownDeath == 1 & Target$NeedsImputation) ||
      length(unique(Target$Horizon)) != 1 || any(table(Target$Treatment) < 2) ||
      !setequal(Target$Treatment, 0:1)) stop("Invalid application imputation inputs.")
  if (Cutpoints[1] != 0 || any(diff(Cutpoints) <= 0) || tail(Cutpoints, 1) < Target$Horizon[1]) {
    stop("Baseline intervals do not cover the imputation horizon.")
  }
  invisible(M)
}


DrawTerminalFollowup <- function(Start, Horizon, Rate, Cutpoints, Uniform) {
  # Invert the conditional cumulative hazard starting at the censoring time.
  N <- length(Start)
  if (!is.matrix(Rate) || nrow(Rate) != N || ncol(Rate) != length(Cutpoints) - 1 ||
      any(!is.finite(Rate)) || any(Rate < 0) || length(Uniform) != N ||
      any(Uniform <= 0 | Uniform >= 1)) stop("Invalid terminal-event sampling inputs.")
  End <- rep(Horizon, N)
  Event <- rep(0, N)
  Remaining <- -log1p(-Uniform)
  for (K in seq_len(ncol(Rate))) {
    Left <- pmax(Start, Cutpoints[K])
    Width <- pmax(0, pmin(Horizon, Cutpoints[K + 1]) - Left)
    Increment <- Rate[, K]*Width
    Hit <- Event == 0 & Width > 0 & Rate[, K] > 0 & Remaining <= Increment
    End[Hit] <- Left[Hit] + Remaining[Hit]/Rate[Hit, K]
    Event[Hit] <- 1
    Alive <- Event == 0
    Remaining[Alive] <- Remaining[Alive] - Increment[Alive]
  }
  if (any(End < Start | End > Horizon)) stop("An imputed death is outside its risk window.")
  list(End = End, Death = Event)
}


GetFutureMoments <- function(Start, Horizon, RecurrenceRate, DeathRate, Cutpoints) {
  # Integrate recurrence intensity times conditional survival, interval by interval.
  Survival <- rep(1, length(Start))
  Count <- Life <- rep(0, length(Start))
  for (K in seq_len(ncol(DeathRate))) {
    Width <- pmax(0, pmin(Horizon, Cutpoints[K + 1]) - pmax(Start, Cutpoints[K]))
    D <- DeathRate[, K]
    AliveTime <- Width
    Positive <- D > 0
    AliveTime[Positive] <- -expm1(-D[Positive]*Width[Positive])/D[Positive]
    AliveTime <- Survival*AliveTime
    Count <- Count + RecurrenceRate[, K]*AliveTime
    Life <- Life + AliveTime
    Survival <- Survival*exp(-D*Width)
  }
  list(Count = Count, RemainingLife = Life, Death = 1 - Survival)
}


ImputeMiBatch <- function(Indices, Draws, Target, Cutpoints, Scenarios, Seed) {
  N <- nrow(Target)
  Need <- which(Target$NeedsImputation)
  M <- length(Indices)
  NewList <- function(Value) setNames(lapply(Scenarios$Scenario, function(Name) matrix(Value, N, M)), Scenarios$Scenario)
  Counts <- ExpectedCounts <- ConditionalCounts <- NewList(Target$ObservedCount)
  Death <- ExpectedDeath <- NewList(Target$KnownDeath)
  End <- ExpectedEnd <- NewList(Target$ObservedEnd)
  DeathRules <- unique(Scenarios$DeathBorrowing)
  if (length(Need)) for (J in seq_along(Indices)) {
    S <- Indices[J]
    # The same independent death/count uniforms couple all scenarios.
    # Draw-index seeds preserve results across worker counts and scenario order.
    SetSimulationSeed(Seed + S)
    UniformDeath <- pmax(.Machine$double.eps, stats::runif(length(Need)))
    UniformCount <- pmax(.Machine$double.eps, stats::runif(length(Need)))
    P <- Target[Need, ]
    Beta <- Draws$Fixed[S, ]
    LogR <- outer(Draws$RecurrenceFrailty[Need, S] +
                    P$Treatment*Beta["Treatment_S1"] + P$Age60*Beta["Age60_S1"],
                  Draws$RecurrenceBaseline[, S], "+")
    LogD <- outer(Draws$DeathFrailty[Need, S] +
                    P$Treatment*Beta["Treatment_S2"] + P$Age60*Beta["Age60_S2"],
                  Draws$DeathBaseline[, S], "+")
    RateR <- exp(LogR)
    AssignedRateD <- exp(LogD)
    if (any(!is.finite(c(RateR, AssignedRateD)))) stop("A sampled future intensity is nonfinite.")
    # Reuse terminal completions within each mortality rule, not across different rules.
    for (Borrowing in DeathRules) {
      RateD <- AssignedRateD*exp(-Borrowing*P$Treatment*Beta["Treatment_S2"])
      if (any(!is.finite(RateD))) stop("A switched death hazard is nonfinite.")
      Terminal <- DrawTerminalFollowup(P$ObservedEnd, P$Horizon[1], RateD, Cutpoints, UniformDeath)
      Exposure <- GetIntervalExposure(P$ObservedEnd, Terminal$End, Cutpoints)
      Moment <- GetFutureMoments(P$ObservedEnd, P$Horizon[1], RateR, RateD, Cutpoints)
      for (I in which(Scenarios$DeathBorrowing == Borrowing)) {
        Name <- Scenarios$Scenario[I]
        Factor <- exp(-Scenarios$RecurrenceBorrowing[I]*P$Treatment*Beta["Treatment_S1"])
        Mean <- rowSums(Exposure*RateR)*Factor
        Missing <- stats::qpois(UniformCount, Mean)
        if (any(!is.finite(c(Mean, Missing))) || any(Mean < 0 | Missing < 0) ||
            any(Missing != floor(Missing)) || any(Missing[rowSums(Exposure) == 0] != 0)) {
          stop("A predictive hospitalization count is invalid.")
        }
        Counts[[Name]][Need, J] <- P$ObservedCount + Missing
        ConditionalCounts[[Name]][Need, J] <- P$ObservedCount + Mean
        ExpectedCounts[[Name]][Need, J] <- P$ObservedCount + Moment$Count*Factor
        Death[[Name]][Need, J] <- Terminal$Death
        End[[Name]][Need, J] <- Terminal$End
        ExpectedDeath[[Name]][Need, J] <- Moment$Death
        ExpectedEnd[[Name]][Need, J] <- P$ObservedEnd + Moment$RemainingLife
      }
    }
  }
  list(Indices = Indices, Counts = Counts, ExpectedCounts = ExpectedCounts,
       ConditionalCounts = ConditionalCounts, Death = Death, End = End,
       ExpectedDeath = ExpectedDeath, ExpectedEnd = ExpectedEnd)
}


GetCompletedAnalyses <- function(Values, Expected, Target, Scenario, Outcome, Meta) {
  Active <- Target$Treatment == 1
  N1 <- sum(Active)
  N0 <- sum(!Active)
  Mean1 <- colMeans(Values[Active, , drop = FALSE])
  Mean0 <- colMeans(Values[!Active, , drop = FALSE])
  V1 <- colSums(sweep(Values[Active, , drop = FALSE], 2, Mean1)^2)/(N1*(N1 - 1))
  V0 <- colSums(sweep(Values[!Active, , drop = FALSE], 2, Mean0)^2)/(N0*(N0 - 1))
  data.frame(Scenario = Scenario, Outcome = Outcome, Imputation = Meta$Imputation, Chain = Meta$Chain, Iteration = Meta$Iteration,
             Estimate = Mean1 - Mean0, WithinVariance = V1 + V0,
             ActiveMean = Mean1, ControlMean = Mean0, ActiveVariance = V1, ControlVariance = V0,
             ActiveParticipants = N1, ControlParticipants = N0, NoMissing = !any(Target$NeedsImputation),
             ExpectedEstimate = colMeans(Expected[Active, , drop = FALSE]) - colMeans(Expected[!Active, , drop = FALSE]))
}


PoolMiScalar <- function(Data, ConfidenceLevel = 0.95) {
  M <- nrow(Data)
  Required <- c("Estimate", "WithinVariance", "ActiveVariance", "ControlVariance", "ExpectedEstimate")
  if (M < 2 || any(!is.finite(as.matrix(Data[, Required]))) ||
      any(Data$ActiveVariance < 0 | Data$ControlVariance < 0) ||
      any(abs(Data$WithinVariance - Data$ActiveVariance - Data$ControlVariance) > 1e-10) ||
      length(unique(Data$Scenario)) != 1 || length(unique(Data$Outcome)) != 1 ||
      anyDuplicated(Data$Imputation)) stop("Invalid completed-data analyses.")
  Q <- mean(Data$Estimate)
  U <- mean(Data$WithinVariance)
  B <- stats::var(Data$Estimate)
  T <- U + (1 + 1/M)*B
  V1 <- mean(Data$ActiveVariance)
  V0 <- mean(Data$ControlVariance)
  N1 <- Data$ActiveParticipants[1]
  N0 <- Data$ControlParticipants[1]
  CompleteDf <- if (U > 0) U^2/(V1^2/(N1 - 1) + V0^2/(N0 - 1)) else NA_real_
  Lambda <- if (T > 0) (1 + 1/M)*B/T else NA_real_
  OldDf <- if (is.finite(Lambda) && Lambda > 0) (M - 1)/Lambda^2 else Inf
  ObservedDf <- if (is.finite(CompleteDf) && is.finite(Lambda)) {
    CompleteDf*(CompleteDf + 1)/(CompleteDf + 3)*(1 - Lambda)
  } else NA_real_
  Df <- if (isTRUE(Data$NoMissing[1])) CompleteDf else 1/(1/OldDf + 1/ObservedDf)
  if (isTRUE(Data$NoMissing[1]) && B > 1e-20) stop("Complete outcomes have nonzero imputation variance.")
  Valid <- is.finite(T) && T > 0 && U > 0 && is.finite(Df) && Df > 0
  Se <- sqrt(T)
  Critical <- if (Valid) stats::qt((1 + ConfidenceLevel)/2, Df) else NA_real_
  Residual <- Data$Estimate - Data$ExpectedEstimate
  data.frame(Scenario = Data$Scenario[1], Outcome = Data$Outcome[1], Imputations = M,
             ActiveMean = mean(Data$ActiveMean), ControlMean = mean(Data$ControlMean), Estimate = Q,
             RubinSe = Se, Lower = Q - Critical*Se, Upper = Q + Critical*Se,
             WithinVariance = U, BetweenVariance = B, TotalVariance = T,
             CompleteDf = CompleteDf, OldDf = OldDf, ObservedDf = ObservedDf, Df = Df,
             Lambda = Lambda, PoolingValid = Valid,
             NumericalMcse = GetChainMcse(Data$Estimate, Data),
             NumericalVarianceFraction = if (T > 0) GetChainMcse(Data$Estimate, Data)^2/T else NA_real_,
             ExpectedEstimate = mean(Data$ExpectedEstimate), ExpectedMcse = GetChainMcse(Data$ExpectedEstimate, Data),
             PredictiveSamplingDifference = mean(Residual), PredictiveSamplingMcse = GetChainMcse(Residual, Data),
             ConfidenceLevel = ConfidenceLevel, NoMissing = Data$NoMissing[1], row.names = NULL)
}


PoolMiAnalyses <- function(Analyses, ConfidenceLevel) {
  Groups <- unique(Analyses[, c("Scenario", "Outcome")])
  do.call(rbind, lapply(seq_len(nrow(Groups)), function(I) {
    Keep <- Analyses$Scenario == Groups$Scenario[I] & Analyses$Outcome == Groups$Outcome[I]
    PoolMiScalar(Analyses[Keep, ], ConfidenceLevel)
  }))
}


CheckPlotSettings <- function(PlotSettings) {
  Default <- HfActionPlotSettings
  if (!is.list(PlotSettings) || any(!names(PlotSettings) %in% names(Default))) stop("Unknown plot setting.")
  P <- utils::modifyList(Default, PlotSettings)
  Numeric <- names(Default)[vapply(Default, is.numeric, logical(1))]
  for (Name in Numeric) if (!is.numeric(P[[Name]]) || any(!is.finite(P[[Name]])) || any(P[[Name]] <= 0)) {
    stop(paste("Invalid plot setting:", Name))
  }
  for (Name in setdiff(Numeric, c("MainMargins", "ModelMargins"))) if (length(P[[Name]]) != 1) stop("A plot size must be a single number.")
  if (length(P$MainMargins) != 4 || sum(P$MainMargins[c(2, 4)]) >= P$FigureWidth ||
      sum(P$MainMargins[c(1, 3)]) >= P$FigureHeight) stop("Plot margins leave no drawing area.")
  ModelPanelHeight <- min((P$ModelObservedFitHeight - P$ModelLegendSpace)/2,
                          P$ModelBaselineHeight, P$ModelCoxSnellHeight - P$ModelLegendSpace)
  if (length(P$ModelMargins) != 4 || sum(P$ModelMargins[c(2, 4)]) >= P$ModelFigureWidth/2 ||
      sum(P$ModelMargins[c(1, 3)]) >= ModelPanelHeight) {
    stop("Model figure margins leave no drawing area.", call. = FALSE)
  }
  if (length(P$FontFamily) != 1 || !nzchar(P$FontFamily) || length(P$PdfFallbackFamily) != 1 ||
      !nzchar(P$PdfFallbackFamily) || length(P$ScenarioColours) < 1 || length(P$ArmColours) != 2) stop("Invalid plot fonts or colours.")
  tryCatch(grDevices::col2rgb(c(P$ScenarioColours, P$ArmColours, P$ReferenceColour)),
           error = function(E) stop("Invalid plot colour."))
  P
}


OpenApplicationPdf <- function(Path, P, Width, Height) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(Path, width = Width, height = Height, pointsize = P$BaseFontSize,
                         family = P$FontFamily, onefile = TRUE)
  } else {
    grDevices::pdf(Path, width = Width, height = Height, pointsize = P$BaseFontSize,
                   family = P$PdfFallbackFamily, onefile = TRUE, useDingbats = FALSE, title = "")
  }
}


SetApplicationPanels <- function(P, Rows = 1, Columns = 1, Main = FALSE) {
  graphics::par(mfrow = c(Rows, Columns), cex = 1, bty = "n", las = 1,
                cex.axis = P$TickLabelSize/P$BaseFontSize, cex.lab = P$AxisLabelSize/P$BaseFontSize,
                mgp = c(2.8, 0.55, 0), tcl = -0.25, oma = c(0, 0, 0, 0),
                mai = if (Main) P$MainMargins else c(1.0, 0.82, 0.45, 0.25))
}


DrawApplicationAxes <- function(P, XLabel = "", YLabel = "", YAt = NULL, YLabels = NULL,
                                XAt = NULL, XLabels = NULL, Panel = NULL) {
  # The surrounding frame is absent. Restore the bottom and left axes explicitly.
  if (is.null(XAt)) graphics::axis(1, lwd = 0, lwd.ticks = P$AxisLineWidth) else {
    graphics::axis(1, at = XAt, labels = if (is.null(XLabels)) XAt else XLabels, lwd = 0, lwd.ticks = P$AxisLineWidth)
  }
  if (is.null(YAt)) graphics::axis(2, lwd = 0, lwd.ticks = P$AxisLineWidth) else {
    graphics::axis(2, at = YAt, labels = if (is.null(YLabels)) YAt else YLabels,
                   cex.axis = P$ScenarioLabelSize/P$BaseFontSize, lwd = 0, lwd.ticks = P$AxisLineWidth)
  }
  Region <- graphics::par("usr")
  if (graphics::par("xlog")) Region[1:2] <- 10^Region[1:2]
  if (graphics::par("ylog")) Region[3:4] <- 10^Region[3:4]
  graphics::segments(Region[1], Region[3], Region[2], Region[3], lwd = P$AxisLineWidth, xpd = NA)
  graphics::segments(Region[1], Region[3], Region[1], Region[4], lwd = P$AxisLineWidth, xpd = NA)
  graphics::mtext(XLabel, side = 1, line = P$XLabelLine, cex = P$AxisLabelSize/P$BaseFontSize)
  graphics::mtext(YLabel, side = 2, line = P$YLabelLine, las = 0, cex = P$AxisLabelSize/P$BaseFontSize)
  if (!is.null(Panel)) graphics::mtext(Panel, side = 3, line = 0.6, adj = 0,
                                       cex = P$PanelLabelSize/P$BaseFontSize, font = 2)
}


DrawApplicationPlot <- function(P, X, Y, XLabel = "", YLabel = "", Panel = NULL, ...) {
  graphics::plot(X, Y, axes = FALSE, ann = FALSE, bty = "n", ...)
  DrawApplicationAxes(P, XLabel, YLabel, Panel = Panel)
}


DrawApplicationLegend <- function(P, Labels, Colours, Types = NULL, Symbols = NULL, Fill = NULL, Columns = length(Labels)) {
  # Figure coordinates keep the legend below the axes for linear and log plots.
  X <- graphics::grconvertX(mean(graphics::par("plt")[1:2]), from = "nfc", to = "user")
  Y <- graphics::grconvertY(graphics::par("plt")[3], from = "nfc", to = "inches") - P$LegendBelowInches
  Y <- graphics::grconvertY(Y, from = "inches", to = "user")
  Arguments <- list(x = X, y = Y, legend = Labels, col = Colours, bty = "n", xpd = NA,
                    xjust = 0.5, yjust = 1, ncol = Columns, cex = P$LegendSize/P$BaseFontSize)
  if (!is.null(Types)) Arguments <- c(Arguments, list(lty = Types, lwd = P$LineWidth))
  if (!is.null(Symbols)) Arguments$pch <- Symbols
  if (!is.null(Fill)) Arguments <- c(Arguments, list(fill = Fill, border = NA))
  do.call(graphics::legend, Arguments)
}


GetScenarioLabels <- function(Scenarios, Wrap = FALSE) {
  Labels <- as.character(Scenarios$Scenario)
  Map <- c(MAR = "Assigned treatment", J2RHospitalization = "Hospitalization switch",
           J2RDeath = "Mortality switch", J2RBoth = "Both switch",
           WJ2R50 = "50% hospitalization borrowing", J2R = "Hospitalization switch")
  Match <- match(Labels, names(Map))
  Labels[!is.na(Match)] <- unname(Map[Match[!is.na(Match)]])
  if (Wrap) Labels <- vapply(Labels, function(Label) paste(strwrap(Label, width = 20), collapse = "\n"), character(1))
  Labels
}


WriteMiPlots <- function(Tables, Settings, MainPath, SupplementPath, PlotSettings) {
  P <- CheckPlotSettings(PlotSettings)
  Names <- as.character(Settings$Scenarios$Scenario)
  Labels <- GetScenarioLabels(Settings$Scenarios)
  Colours <- rep(P$ScenarioColours, length.out = length(Names))
  Paths <- c(HospitalizationCount = MainPath,
             DeathRisk = file.path(dirname(MainPath), "MortalitySensitivity.pdf"),
             RestrictedMeanSurvival = file.path(dirname(MainPath), "SurvivalSensitivity.pdf"))
  Scales <- c(HospitalizationCount = 1, DeathRisk = 100, RestrictedMeanSurvival = 365.25)
  AxisLabels <- c(HospitalizationCount = "Difference in mean hospitalizations",
                  DeathRisk = "Difference in death risk (percentage points)",
                  RestrictedMeanSurvival = "Difference in restricted mean survival (days)")
  Device <- NULL
  on.exit(if (!is.null(Device) && Device %in% grDevices::dev.list()) grDevices::dev.off(Device), add = TRUE)
  for (Outcome in names(Paths)) {
    X <- Tables$MiPooled[Tables$MiPooled$Outcome == Outcome, ]
    X <- X[match(Names, X$Scenario), ]
    if (anyNA(X$Scenario)) stop("A sensitivity scenario is missing from the pooled results.")
    for (Field in c("Lower", "Upper", "Estimate")) X[[Field]] <- X[[Field]]*Scales[Outcome]
    OpenApplicationPdf(Paths[Outcome], P, P$FigureWidth, max(P$FigureHeight, 1.8 + 0.45*nrow(X)))
    Device <- grDevices::dev.cur()
    SetApplicationPanels(P, Main = TRUE)
    Values <- c(X$Lower, X$Upper, X$Estimate, 0)
    Limits <- range(Values[is.finite(Values)])
    Padding <- max(diff(Limits)*0.08, 0.01)
    Y <- rev(seq_len(nrow(X)))
    graphics::plot(X$Estimate, Y, type = "n", axes = FALSE, ann = FALSE,
                   xlim = Limits + c(-Padding, Padding), ylim = c(0.5, nrow(X) + 0.5), xaxs = "i", yaxs = "i", bty = "n")
    graphics::abline(v = 0, lty = 2, col = P$ReferenceColour, lwd = P$ReferenceLineWidth)
    graphics::segments(X$Lower, Y, X$Upper, Y, col = Colours, lwd = P$LineWidth)
    graphics::segments(X$Lower, Y - 0.055, X$Lower, Y + 0.055, col = Colours, lwd = P$LineWidth)
    graphics::segments(X$Upper, Y - 0.055, X$Upper, Y + 0.055, col = Colours, lwd = P$LineWidth)
    graphics::points(X$Estimate, Y, pch = 16, col = Colours, cex = P$PointSize)
    DrawApplicationAxes(P, paste0(AxisLabels[Outcome], "\nExercise minus usual care"), YAt = Y, YLabels = Labels)
    grDevices::dev.off(Device)
  }
  OpenApplicationPdf(SupplementPath, P, P$SupplementSingleWidth, P$SupplementSingleHeight)
  Device <- grDevices::dev.cur()
  BurdenFields <- c("MeanImputedHospitalizations", "MeanAdditionalDeaths", "MeanAdditionalPersonYears")
  BurdenLabels <- c("Mean imputed hospitalizations per participant", "Mean imputed deaths per participant",
                    "Mean imputed time alive per participant (days)")
  for (I in seq_along(BurdenFields)) {
    SetApplicationPanels(P)
    Heights <- vapply(Names, function(Name) {
      Z <- Tables$MiBurden[Tables$MiBurden$Scenario == Name, ]
      Z[[BurdenFields[I]]][match(0:1, Z$Treatment)]/Z$Participants[match(0:1, Z$Treatment)]
    }, numeric(2))
    if (I == 3) Heights <- Heights*365.25
    At <- graphics::barplot(Heights, beside = TRUE, col = P$ArmColours, border = NA,
                            axes = FALSE, axisnames = FALSE, ylim = c(0, max(Heights, 0.01)*1.12), yaxs = "i")
    DrawApplicationAxes(P, YLabel = BurdenLabels[I], XAt = colMeans(At), XLabels = FALSE)
    LabelY <- graphics::grconvertY(graphics::par("plt")[3], from = "nfc", to = "inches") - P$ScenarioLabelOffsetInches
    LabelY <- graphics::grconvertY(LabelY, from = "inches", to = "user")
    graphics::text(colMeans(At), LabelY, GetScenarioLabels(Settings$Scenarios, Wrap = TRUE),
                   adj = c(0.5, 1), cex = P$ScenarioLabelSize/P$BaseFontSize, xpd = NA)
    DrawApplicationLegend(P, c("Usual care", "Exercise"), P$ArmColours, Fill = P$ArmColours)
  }
  for (Outcome in names(Paths)) {
    SetApplicationPanels(P)
    graphics::par(mai = c(1.35, 0.9, 0.35, 0.25))
    Z <- Tables$MiPrecision[Tables$MiPrecision$Outcome == Outcome, ]
    Z$Estimate <- Z$Estimate*Scales[Outcome]
    DrawApplicationPlot(P, range(Z$Imputations), range(Z$Estimate), "Number of imputations",
                        AxisLabels[Outcome], type = "n")
    for (I in seq_along(Names)) {
      A <- Z[Z$Scenario == Names[I], ]
      A <- A[order(A$Imputations), ]
      graphics::lines(A$Imputations, A$Estimate, type = "b", pch = 16, col = Colours[I],
                      lwd = P$LineWidth, cex = P$PointSize)
    }
    DrawApplicationLegend(P, Labels, Colours, Types = rep(1, length(Names)), Symbols = rep(16, length(Names)), Columns = 2)
  }
  grDevices::dev.off(Device)
  invisible(c(Paths, SupplementPath))
}


WriteManuscriptModelPlots <- function(Tables, Directory, PlotSettings = HfActionPlotSettings) {
  P <- CheckPlotSettings(PlotSettings)
  Required <- list(
    RiskSetCalibration = c("Endpoint", "Treatment", "Start", "Stop", "Observed", "Mean", "Lower", "Upper"),
    BaselineHazards = c("Endpoint", "Start", "Stop", "Median", "Lower", "Upper"),
    CoxSnellCurves = c("Endpoint", "Group", "Source", "Residual", "AtRisk", "CumulativeHazard")
  )
  for (Name in names(Required)) {
    if (!is.data.frame(Tables[[Name]]) || !all(Required[[Name]] %in% names(Tables[[Name]]))) {
      stop(paste("Missing model plot data:", Name), call. = FALSE)
    }
  }
  if (!dir.exists(Directory) && !dir.create(Directory, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create the manuscript plot folder.", call. = FALSE)
  }
  
  Q <- P
  Q$AxisLabelSize <- P$ModelAxisLabelSize
  Q$TickLabelSize <- P$ModelTickLabelSize
  Q$LegendSize <- P$ModelLegendSize
  Q$XLabelLine <- P$ModelXLabelLine
  Q$YLabelLine <- P$ModelYLabelLine
  Q$LineWidth <- P$ModelLineWidth
  Q$PointSize <- P$ModelPointSize
  Paths <- file.path(Directory, c("ModelObservedFit.pdf", "ModelBaselineHazards.pdf", "ModelCoxSnell.pdf"))
  Names <- c("Recurrence", "Death")
  Colours <- P$ArmColours
  
  SetPanels <- function(Rows, Columns, Legend = FALSE) {
    graphics::par(mfrow = c(Rows, Columns), cex = 1, bty = "n", las = 1,
                  cex.axis = Q$TickLabelSize/Q$BaseFontSize,
                  cex.lab = Q$AxisLabelSize/Q$BaseFontSize, mgp = c(2.8, 0.55, 0),
                  tcl = -0.25, xpd = FALSE, oma = c(0, 0, 0, 0),
                  omi = c(if (Legend) P$ModelLegendSpace else 0, 0, 0, 0),
                  mai = P$ModelMargins)
  }
  
  DrawLegend <- function(Height, Labels, LineColours, Types, Symbols = NULL) {
    # Device coordinates place one shared legend below every panel.
    X <- graphics::grconvertX(0.5, from = "ndc", to = "user")
    Y <- graphics::grconvertY(P$ModelLegendSpace/(2*Height), from = "ndc", to = "user")
    Arguments <- list(x = X, y = Y, legend = Labels, col = LineColours,
                      lty = Types, lwd = Q$LineWidth, bty = "n", xpd = NA,
                      xjust = 0.5, yjust = 0.5, ncol = length(Labels),
                      cex = Q$LegendSize/Q$BaseFontSize)
    if (!is.null(Symbols)) Arguments$pch <- Symbols
    do.call(graphics::legend, Arguments)
  }
  
  DrawObserved <- function() {
    OpenApplicationPdf(Paths[1], Q, P$ModelFigureWidth, P$ModelObservedFitHeight)
    on.exit(grDevices::dev.off(), add = TRUE)
    SetPanels(2, 2, Legend = TRUE)
    for (Endpoint in Names) for (Arm in 0:1) {
      D <- Tables$RiskSetCalibration
      D <- D[D$Endpoint == Endpoint & D$Treatment == Arm, ]
      D <- D[order(D$Start, D$Stop), ]
      Values <- unlist(D[c("Start", "Stop", "Observed", "Mean", "Lower", "Upper")], use.names = FALSE)
      if (!nrow(D) || any(!is.finite(Values))) stop("Observed fit summaries are incomplete.", call. = FALSE)
      Mid <- (D$Start + D$Stop)/2
      ArmLabel <- if (Arm == 0) "usual care" else "exercise"
      YLabel <- paste0(if (Endpoint == "Recurrence") "Hospitalizations" else "Deaths", " (", ArmLabel, ")")
      DrawApplicationPlot(Q, Mid, D$Observed, "Follow-up (years)", YLabel,
                          pch = 16, col = "black", cex = Q$PointSize,
                          xlim = range(c(D$Start, D$Stop)), ylim = c(0, max(c(D$Observed, D$Upper))*1.08 + 1))
      graphics::segments(Mid, D$Lower, Mid, D$Upper, col = Colours[1], lwd = Q$LineWidth)
      graphics::lines(Mid, D$Mean, col = Colours[1], type = "b", pch = 1,
                      cex = Q$PointSize, lwd = Q$LineWidth)
    }
    DrawLegend(P$ModelObservedFitHeight, c("Observed", "Fitted intensity"),
               c("black", Colours[1]), c(0, 1), c(16, 1))
  }
  
  DrawBaselines <- function() {
    OpenApplicationPdf(Paths[2], Q, P$ModelFigureWidth, P$ModelBaselineHeight)
    on.exit(grDevices::dev.off(), add = TRUE)
    SetPanels(1, 2)
    for (Endpoint in Names) {
      D <- Tables$BaselineHazards
      D <- D[D$Endpoint == Endpoint, ]
      D <- D[order(D$Start, D$Stop), ]
      Values <- unlist(D[c("Start", "Stop", "Median", "Lower", "Upper")], use.names = FALSE)
      if (!nrow(D) || any(!is.finite(Values))) stop("Baseline summaries are incomplete.", call. = FALSE)
      Mid <- (D$Start + D$Stop)/2
      YLabel <- if (Endpoint == "Recurrence") "Baseline recurrence intensity (per year)" else "Baseline death hazard (per year)"
      DrawApplicationPlot(Q, Mid, D$Median, "Follow-up (years)", YLabel,
                          type = "b", pch = 16, col = Colours[1], lwd = Q$LineWidth, cex = Q$PointSize,
                          xlim = range(c(D$Start, D$Stop)), ylim = c(0, max(D$Upper)*1.06))
      graphics::segments(Mid, D$Lower, Mid, D$Upper, col = Colours[1], lwd = Q$LineWidth)
    }
  }
  
  DrawCoxSnell <- function() {
    OpenApplicationPdf(Paths[3], Q, P$ModelFigureWidth, P$ModelCoxSnellHeight)
    on.exit(grDevices::dev.off(), add = TRUE)
    SetPanels(1, 2, Legend = TRUE)
    Sources <- "Mcmc"
    for (Endpoint in Names) {
      D <- Tables$CoxSnellCurves
      D <- D[D$Endpoint == Endpoint & D$Group == "All" & D$Source %in% Sources, ]
      if (!nrow(D) || !all(Sources %in% D$Source) ||
          any(!is.finite(unlist(D[c("Residual", "AtRisk", "CumulativeHazard")], use.names = FALSE)))) {
        stop("Cox-Snell summaries are incomplete.", call. = FALSE)
      }
      Supported <- vapply(Sources, function(Source) {
        At <- D$Residual[D$Source == Source & D$AtRisk >= 10]
        if (length(At)) max(At) else 0.1
      }, numeric(1))
      Limit <- max(0.1, min(Supported))
      YLimit <- max(Limit, D$CumulativeHazard[D$Residual <= Limit])
      XLabel <- if (Endpoint == "Recurrence") "Recurrence cumulative intensity" else "Death cumulative hazard"
      DrawApplicationPlot(Q, 0, 0, XLabel, "Nelson-Aalen cumulative hazard",
                          type = "n", xlim = c(0, Limit), ylim = c(0, YLimit))
      graphics::abline(0, 1, lty = 3, col = P$ReferenceColour, lwd = P$ReferenceLineWidth)
      for (Index in seq_along(Sources)) {
        Curve <- D[D$Source == Sources[Index], ]
        Curve <- Curve[order(Curve$Residual), ]
        graphics::lines(Curve$Residual, Curve$CumulativeHazard, type = "s",
                        col = Colours[Index], lty = Index, lwd = Q$LineWidth)
      }
    }
    DrawLegend(P$ModelCoxSnellHeight, c("Posterior mean rates", "Exponential reference"), c(Colours[1], P$ReferenceColour), c(1, 3))
  }
  
  DrawObserved()
  DrawBaselines()
  DrawCoxSnell()
  invisible(Paths)
}


LatexEscape <- function(Text) {
  Map <- c("\\" = "\\textbackslash{}", "&" = "\\&", "%" = "\\%", "$" = "\\$",
           "#" = "\\#", "_" = "\\_", "{" = "\\{", "}" = "\\}", "~" = "\\textasciitilde{}", "^" = "\\textasciicircum{}")
  vapply(strsplit(as.character(Text), "", fixed = TRUE), function(Characters) {
    At <- match(Characters, names(Map))
    Characters[!is.na(At)] <- Map[At[!is.na(At)]]
    paste(Characters, collapse = "")
  }, character(1), USE.NAMES = FALSE)
}


LatexNumber <- function(Value, Digits = 3) {
  Value[is.finite(Value) & abs(Value) < 0.5*10^(-Digits)] <- 0
  ifelse(is.finite(Value), formatC(Value, format = "f", digits = Digits, decimal.mark = "."), "{NA}")
}


LatexColumn <- function(Values, Digits) {
  Values <- Values[is.finite(Values)]
  if (!length(Values)) Values <- 0
  Width <- max(1, nchar(as.character(floor(abs(round(Values, Digits))))))
  paste0("S[table-format=", if (any(Values < 0)) "-" else "", Width, ".", Digits, "]")
}


LatexBound <- function(Values, Digits, Prefix = "", Suffix = "") {
  vapply(Values, function(Value) {
    if (!is.finite(Value)) return(paste0("\\mbox{", Prefix, "\\textnormal{NA}", Suffix, "}"))
    Number <- LatexNumber(Value, Digits)
    Column <- LatexColumn(as.numeric(Number), Digits)
    Specification <- sub("^S\\[table-format=(.*)\\]$", "\\1", Column)
    paste0("\\mbox{", Prefix, "\\tablenum[table-format=", Specification, "]{", Number, "}", Suffix, "}")
  }, character(1), USE.NAMES = FALSE)
}


WriteLatexTable <- function(Path, Caption, Label, Columns, Header, Rows, Footnote, TableSettings) {
  FontCommands <- c("tiny", "scriptsize", "footnotesize", "small", "normalsize", "large", "Large")
  if (!all(c(TableSettings$TableFont, TableSettings$CaptionFont, TableSettings$FootnoteFont) %in% FontCommands) ||
      any(!is.finite(c(TableSettings$ColumnSepPt, TableSettings$RowStretch))) ||
      any(c(TableSettings$ColumnSepPt, TableSettings$RowStretch) <= 0)) stop("Invalid LaTeX table settings.")
  Lines <- c("% Generated from the saved pooled results. Requires LatexPreamble.tex.",
             "\\begin{table}[!htbp]", "\\begingroup", "\\centering",
             paste0("\\", TableSettings$TableFont),
             paste0("\\setlength{\\tabcolsep}{", TableSettings$ColumnSepPt, "pt}"),
             paste0("\\renewcommand{\\arraystretch}{", TableSettings$RowStretch, "}"),
             paste0("\\captionsetup{justification=justified,singlelinecheck=false,font=", TableSettings$CaptionFont, "}"),
             paste0("\\caption{", Caption, "}"), paste0("\\label{", Label, "}"),
             paste0("\\begin{tabular}{@{}", Columns, "@{}}"), "\\toprule", Header, "\\midrule", Rows,
             "\\bottomrule", "\\end{tabular}", "\\par\\vspace{4pt}", "\\begin{minipage}{\\linewidth}",
             paste0("\\", TableSettings$FootnoteFont, "\\justifying"), Footnote,
             "\\end{minipage}", "\\endgroup", "\\end{table}")
  writeLines(Lines, Path, useBytes = TRUE)
  invisible(Path)
}


WriteApplicationTables <- function(Tables, Settings, Directory, TableSettings) {
  Digits <- c(TableSettings$CountDigits, TableSettings$SurvivalDigits)
  if (length(Digits) != 2 || any(!is.finite(Digits)) || any(Digits < 0 | Digits > 8 | Digits != floor(Digits))) {
    stop("Use 0 to 8 decimal places in the tables.")
  }
  Settings$Scenarios <- CheckMiScenarios(Settings$Scenarios)
  Main <- file.path(Directory, "Manuscript")
  Supplement <- file.path(Directory, "Supplement")
  Row <- function(Cells) paste0(paste(Cells, collapse = " & "), " \\\\")
  H <- Tables$HorizonSummary[match(0:1, Tables$HorizonSummary$Treatment), ]
  Horizon <- format(Settings$Horizon, trim = TRUE)
  Confidence <- format(100*Settings$ConfidenceLevel, trim = TRUE)
  Scenarios <- as.character(Settings$Scenarios$Scenario)
  Labels <- LatexEscape(GetScenarioLabels(Settings$Scenarios))
  Pooled <- Tables$MiPooled
  Counts <- Pooled[Pooled$Outcome == "HospitalizationCount", ]
  Counts <- Counts[match(Scenarios, Counts$Scenario), ]
  ArmColumn <- paste0(LatexColumn(c(H$Participants, H$ObservedHospitalizationsByHorizon, H$ObservedPersonYears), 1),
    "@{\\,}r")
  ObservedRow <- function(Label, Values, Digits = 0, Percent = NULL) {
    Cells <- c(Label)
    for (I in 1:2) Cells <- c(Cells, LatexNumber(Values[I], Digits),
      if (is.null(Percent)) "" else LatexBound(Percent[I], 1, Prefix = "(", Suffix = ")"))
    Row(Cells)
  }
  Rows <- c(ObservedRow("Participants", H$Participants),
            ObservedRow("Observed hospitalizations", H$ObservedHospitalizationsByHorizon),
            ObservedRow("Observed deaths", H$ObservedDeathsByHorizon),
            ObservedRow("Observed person-years", H$ObservedPersonYears, 1),
            ObservedRow("Alive and observed through the horizon", H$ObservedThroughHorizon),
            ObservedRow("Censored alive before the horizon, $n$ (\\%)", H$NeedImputation, Percent = H$PercentNeedImputation))
  WriteLatexTable(file.path(Main, "TableObservedData.tex"),
                  paste0("Observed data through ", Horizon, " years in the public HF-ACTION subset."), "tab:HfObserved",
                  paste0("l", ArmColumn, "@{\\hspace{16pt}}", ArmColumn),
                  " & \\multicolumn{2}{l}{\\textbf{Usual care}} & \\multicolumn{2}{l@{}}{\\textbf{Exercise}} \\\\", Rows,
                  paste0("Only participants censored alive before ", Horizon, " years require imputation. Percentages use the number of participants in each arm as the denominator. ",
                         "Observed hospitalizations and person-time stop at the earliest of recorded censoring, death and the analysis horizon. All available observed follow-up, including follow-up beyond the horizon, enters model fitting. ",
                         "The public dataset comprises participants with nonischemic heart failure and baseline cardiopulmonary exercise test duration of at most 9 minutes; it is a subset of the trial."), TableSettings)
  ResultColumns <- function(X, Digits) paste("l", LatexColumn(X$ControlMean, Digits), LatexColumn(X$ActiveMean, Digits),
    LatexColumn(X$Estimate, Digits), LatexColumn(X$RubinSe, Digits), "r@{,\\,}r")
  ResultRows <- function(X, Labels, Digits) {
    Lower <- LatexBound(X$Lower, Digits, Prefix = "(")
    Upper <- LatexBound(X$Upper, Digits, Suffix = ")")
    vapply(seq_len(nrow(X)), function(I) {
      Row(c(Labels[I], LatexNumber(unlist(X[I, c("ControlMean", "ActiveMean", "Estimate", "RubinSe")]), Digits),
        Lower[I], Upper[I]))
    }, character(1))
  }
  Header <- Row(c("", "\\multicolumn{1}{l}{\\textbf{Usual care}}",
    "\\multicolumn{1}{l}{\\textbf{Exercise}}", "\\multicolumn{1}{l}{\\textbf{Difference}}",
    "\\multicolumn{1}{l}{\\textbf{SE}}", paste0("\\multicolumn{2}{l@{}}{\\textbf{", Confidence, "\\% CI}}")))
  WriteLatexTable(file.path(Main, "TableHospitalizations.tex"),
                  paste0("Mean hospitalization counts through ", Horizon, " years under post-censoring sensitivity assumptions."), "tab:HfHospitalizations",
                  ResultColumns(Counts, TableSettings$CountDigits), Header, ResultRows(Counts, Labels, TableSettings$CountDigits),
                  paste0("Differences are exercise minus usual care. Counts accumulate before death or the horizon, with all participants retained in the denominator. ",
                         "Assigned treatment retains both treatment components after censoring. Other scenarios remove the hospitalization component, the mortality component, or both; see Table~\\ref{tab:HfScenarios}. ",
                         "The participant's covariates and shared frailty are retained. Results use ", Settings$Imputations,
                         " imputations per scenario. SE and CI denote the Rubin standard error and interval, using Welch/Barnard--Rubin degrees of freedom. ",
                         "These are hypothetical assumptions after recorded censoring, not analyses of documented treatment discontinuation. Numerical precision is reported separately from repeated sampling coverage."), TableSettings)
  Rows <- vapply(seq_along(Scenarios), function(I) Row(c(Labels[I],
                                                         LatexNumber(Settings$Scenarios$RecurrenceBorrowing[I], 2),
                                                         LatexNumber(Settings$Scenarios$DeathBorrowing[I], 2))), character(1))
  WriteLatexTable(file.path(Main, "TableScenarios.tex"), "Borrowing assumptions after recorded censoring.", "tab:HfScenarios",
                  "l S[table-format=1.2] S[table-format=1.2]",
                  Row(c("Assumption", "{Recurrence weight}", "{Mortality weight}")), Rows,
                  "A weight of zero retains the assigned-treatment component of the conditional intensity or hazard; a weight of one removes that component in the active arm. Intermediate values remove the specified fraction of the log treatment effect. Covariates and participant frailties are retained, and reference-arm predictions are unchanged. The assigned-treatment scenario is denoted MAR in the machine-readable outputs.", TableSettings)
  Outcomes <- c("DeathRisk", "RestrictedMeanSurvival")
  Survival <- do.call(rbind, lapply(Outcomes, function(Outcome) {
    X <- Pooled[Pooled$Outcome == Outcome, ]
    X[match(Scenarios, X$Scenario), ]
  }))
  Fields <- c("ControlMean", "ActiveMean", "Estimate", "RubinSe", "Lower", "Upper")
  for (Outcome in Outcomes) for (Weight in unique(Settings$Scenarios$DeathBorrowing)) {
    Peers <- Settings$Scenarios$Scenario[Settings$Scenarios$DeathBorrowing == Weight]
    Z <- Pooled[Pooled$Outcome == Outcome & Pooled$Scenario %in% Peers, Fields]
    if (nrow(Z) != length(Peers) || any(vapply(Z, function(V) {
      all(is.na(V)) || (all(is.finite(V)) && max(abs(V - V[1])) <= 1e-10)
    }, logical(1)) == FALSE)) stop("Survival summaries disagree within a mortality rule.")
  }
  Scale <- ifelse(Survival$Outcome == "DeathRisk", 100, 365.25)
  for (Field in Fields) Survival[[Field]] <- Survival[[Field]]*Scale
  Rows <- c("\\multicolumn{7}{@{}l}{\\textit{Death risk (\\%)}} \\\\",
            ResultRows(Survival[Survival$Outcome == "DeathRisk", ], Labels, TableSettings$SurvivalDigits),
            "\\addlinespace", "\\multicolumn{7}{@{}l}{\\textit{Restricted mean survival (days)}} \\\\",
            ResultRows(Survival[Survival$Outcome == "RestrictedMeanSurvival", ], Labels, TableSettings$SurvivalDigits))
  WriteLatexTable(file.path(Main, "TableSurvival.tex"),
                  paste0("Survival outcomes through ", Horizon, " years under each borrowing assumption."), "tab:HfSurvival",
                  ResultColumns(Survival, TableSettings$SurvivalDigits), Header, Rows,
                  paste0("Differences are exercise minus usual care. Death-risk differences are percentage points; restricted mean survival differences are days, using 365.25 days per year. ",
                         "Survival completions are identical only for scenarios with the same mortality borrowing weight. Hospitalization borrowing does not alter the death hazard conditional on the shared frailty. ",
                         "Shorter survival can reduce hospitalization counts, so the two outcomes should be interpreted together. Intervals use the same Rubin procedure as Table~\\ref{tab:HfHospitalizations}."), TableSettings)
  B <- Tables$MiBurden
  B <- B[order(match(B$Scenario, Scenarios), B$Treatment), ]
  Rows <- vapply(seq_len(nrow(B)), function(I) {
    Row(c(Labels[match(B$Scenario[I], Scenarios)], if (B$Treatment[I] == 0) "Usual care" else "Exercise",
          LatexNumber(B$NeedImputation[I], 0), LatexNumber(B$MeanImputedHospitalizations[I], 2),
          LatexNumber(100*B$ImputedFractionOfMeanTotal[I], 1), LatexNumber(B$MeanAdditionalDeaths[I], 2)))
  }, character(1))
  WriteLatexTable(file.path(Supplement, "TableImputationBurden.tex"),
                  paste0("Amount imputed after censoring through ", Horizon, " years."), "tab:HfBurden",
                  paste("ll", LatexColumn(B$NeedImputation, 0), LatexColumn(B$MeanImputedHospitalizations, 2),
                        LatexColumn(100*B$ImputedFractionOfMeanTotal, 1), LatexColumn(B$MeanAdditionalDeaths, 2)),
                  c(Row(c("", "", "{Participants}", "{Hospitalizations}", "{Imputed}", "{Deaths}")),
                    Row(c("Assumption", "Arm", "{requiring MI}", "{added}", "{fraction (\\%)}", "{added}"))), Rows,
                  "Added hospitalizations and deaths are totals averaged over imputations. The imputed hospitalization fraction is the mean number added divided by the mean completed total in that arm. Death completions are shared only within the same mortality borrowing rule. MI denotes multiple imputation.", TableSettings)
  invisible(NULL)
}


ValidateApplicationSettings <- function(S) {
  for (Name in c("Intervals", "Seed", "Chains", "Warmup", "Sampling", "MaximumWarmupBlocks",
                 "MaximumSamplingBlocks", "ImputationBatchSize")) {
    X <- S[[Name]]
    if (length(X) != 1 || !is.numeric(X) || !is.finite(X) || X < 1 || X != floor(X)) stop(paste("Invalid setting:", Name))
  }
  if (S$Chains < 2 || S$Chains > 16 || S$Warmup < 100 || S$Sampling < 100) stop("Use at least two chains and 100 iterations per block.")
  for (Name in c("FixedEffectSd", "AlphaSd", "FrailtyShape", "FrailtyRate", "BaselineLogSd", "MinimumEss", "Horizon")) {
    X <- S[[Name]]
    if (length(X) != 1 || !is.finite(X) || X <= 0) stop(paste("Invalid setting:", Name))
  }
  if (length(S$BaselineLogMean) != 1 || !is.finite(S$BaselineLogMean) ||
      length(S$MaximumRhat) != 1 || !is.finite(S$MaximumRhat) || S$MaximumRhat <= 1 ||
      length(S$CoreFraction) != 1 || !is.finite(S$CoreFraction) || S$CoreFraction <= 0 || S$CoreFraction > 1 ||
      length(S$ConfidenceLevel) != 1 || !is.finite(S$ConfidenceLevel) || S$ConfidenceLevel <= 0 || S$ConfidenceLevel >= 1) stop("Invalid diagnostic or interval settings.")
  for (Name in c("BlockMoves", "InstallMissingPackages")) {
    if (!is.logical(S[[Name]]) || length(S[[Name]]) != 1 || is.na(S[[Name]])) stop(paste("Invalid setting:", Name))
  }
  MaximumDraws <- S$Chains*S$Sampling*S$MaximumSamplingBlocks
  if (S$Seed + MaximumDraws + 100000 >= .Machine$integer.max ||
      S$Warmup + S$Sampling*(S$MaximumWarmupBlocks + S$MaximumSamplingBlocks) >= .Machine$integer.max) stop("Iteration or seed limit is too large.")
  invisible(S)
}

ReportApplicationBlock <- function(B) {
  message("Block ", B$Block, ": ", B$Stage, "; ", B$AssessedPerChain,
    " draws assessed per chain; max Rhat ", format(round(B$MaximumRhat, 4), nsmall = 4),
    "; min bulk/tail ESS ", round(B$MinimumBulkEss), "/", round(B$MinimumTailEss), ".")
}

FitApplication <- function(Prepared, Cells, Cutpoints, Settings, Directory) {
  Specification <- list(Cutpoints = Cutpoints, BaselineTimes = tail(Cutpoints, -1))
  Input <- list(Patients = Prepared$Patients, Cells = Cells)
  D <- BuildCustomData(Input, Specification, Settings)
  Initial <- GetCustomInitial(Input, D)
  Cache <- file.path(dirname(Directory), "SamplerCache")
  CompileCustomSampler(Cache)
  AuditInitial <- JitterCustomInitial(Initial, Settings$Seed + 999, 1)
  Audit <- AuditCustomKernels(D, AuditInitial)
  utils::write.csv(Audit, file.path(Directory, "ConditionalDensityAudit.csv"), row.names = FALSE)
  ParallelChains <- min(Settings$Chains, Settings$CoreBudget)
  State <- new.env(parent = emptyenv())
  State$Cluster <- NULL
  HadState <- exists("CustomWorkerState", envir = .GlobalEnv, inherits = FALSE)
  PreviousState <- if (HadState) get("CustomWorkerState", envir = .GlobalEnv) else NULL
  assign("CustomWorkerState", State, envir = .GlobalEnv)
  on.exit({
    if (!is.null(State$Cluster)) parallel::stopCluster(State$Cluster)
    if (HadState) assign("CustomWorkerState", PreviousState, envir = .GlobalEnv)
    else rm("CustomWorkerState", envir = .GlobalEnv)
  }, add = TRUE)
  if (ParallelChains > 1) {
    State$Cluster <- parallel::makePSOCKcluster(ParallelChains)
    parallel::clusterExport(State$Cluster, c("CustomSamplerCode", "CompileCustomSampler", "SetSimulationSeed", "RunCustomChain"), envir = environment(FitApplication))
    parallel::clusterCall(State$Cluster, function(Path) {
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
      CompileCustomSampler(Path)
      NULL
    }, Cache)
  }
  Jobs <- lapply(seq_len(Settings$Chains), function(J) list(Chain = J, Data = D,
    Initial = JitterCustomInitial(Initial, Settings$Seed + 1000 + J, J),
    Warmup = Settings$Warmup, Sampling = Settings$Sampling, Seed = Settings$Seed + 2000 + J,
    BlockMoves = Settings$BlockMoves, Continuation = NULL))
  message("Sampling ", Settings$Chains, " chains on ", ParallelChains, " core(s); budget ", Settings$CoreBudget,
    "/", Settings$LogicalCores, ". Each block has ", Settings$Sampling, " draws per chain.")
  RecordedWarnings <- character()
  Fit <- withCallingHandlers(RunCustomAutoFit(Jobs, Settings), warning = function(W) {
    RecordedWarnings <<- c(RecordedWarnings, conditionMessage(W))
    invokeRestart("muffleWarning")
  })
  Fit$Warnings <- unique(RecordedWarnings)
  Fit$ReconstructionDifference <- AuditReturnedStates(Fit$Draws, D)
  Fit$DensityAudit <- Audit
  Fit$Data <- D
  Fit$ConditionalKernelHash <- SimulationHash(CustomSamplerCode)
  message("Sampling finished: ", Fit$StopReason, "; ", Fit$DrawsPerChain*Settings$Chains, " retained joint draws.")
  Fit
}

ExtractApplicationDraws <- function(Fit, Target, Cutpoints) {
  A <- Fit$Draws
  Iterations <- dim(A)[1]
  Chains <- dim(A)[2]
  X <- matrix(A, nrow = Iterations*Chains, dimnames = list(NULL, dimnames(A)[[3]]))
  Fixed <- X[, c("BetaR", "GammaR", "BetaD", "GammaD"), drop = FALSE]
  colnames(Fixed) <- c("Treatment_S1", "Age60_S1", "Treatment_S2", "Age60_S2")
  B <- t(X[, paste0("B[", seq_len(nrow(Target)), "]"), drop = FALSE])
  K <- length(Cutpoints) - 1
  list(Fixed = Fixed, RecurrenceBaseline = t(X[, paste0("RecurrenceLogRate[", seq_len(K), "]"), drop = FALSE]),
    DeathBaseline = t(X[, paste0("DeathLogRate[", seq_len(K), "]"), drop = FALSE]),
    RecurrenceFrailty = B, DeathFrailty = sweep(B, 2, X[, "Alpha"], "*"),
    Meta = data.frame(Imputation = seq_len(nrow(X)), Chain = rep(seq_len(Chains), each = Iterations),
      Iteration = rep(seq_len(Iterations), Chains)))
}

GetChainMcse <- function(Values, Meta) {
  if (length(Values) != nrow(Meta) || any(!is.finite(Values))) return(NA_real_)
  if (length(unique(Values)) == 1) return(0)
  ChainSizes <- table(Meta$Chain)
  if (length(ChainSizes) < 2 || length(unique(ChainSizes)) != 1 || min(ChainSizes) < 4 ||
      anyDuplicated(Meta[, c("Chain", "Iteration")])) stop("Monte Carlo errors require equal lengths and unique iterations in each chain.")
  Ordered <- order(Meta$Chain, Meta$Iteration)
  X <- matrix(Values[Ordered], nrow = ChainSizes[1])
  as.numeric(posterior::mcse_mean(X))
}

SummarizeMiBatch <- function(Indices, Draws, Target, Cutpoints, Scenarios, Seed) {
  Z <- ImputeMiBatch(Indices, Draws, Target, Cutpoints, Scenarios, Seed)
  Need <- Target$NeedsImputation
  Known <- !Need
  Control <- Target$Treatment == 0
  Meta <- Draws$Meta[Indices, ]
  Analyses <- Burden <- Participants <- list()
  for (Name in Scenarios$Scenario) {
    C <- Z$Counts[[Name]]
    D <- Z$Death[[Name]]
    T <- Z$End[[Name]]
    Missing <- C - Target$ObservedCount
    if (any(!is.finite(c(C, D, T))) || any(C != floor(C)) || any(Missing < 0) ||
        any(Missing[Known, , drop = FALSE] != 0) || any(!D %in% 0:1) ||
        any(D[Known, , drop = FALSE] != Target$KnownDeath[Known]) ||
        any(T[Known, , drop = FALSE] != Target$ObservedEnd[Known]) ||
        any(T < Target$ObservedEnd | T > Target$Horizon) ||
        any(Missing[T == Target$ObservedEnd] != 0) ||
        any(D[Need, , drop = FALSE] == 0 & T[Need, , drop = FALSE] != Target$Horizon[Need])) stop("An imputation violated an observed history or endpoint.")
    if (!identical(C[Control, , drop = FALSE], Z$Counts$MAR[Control, , drop = FALSE]) ||
        !identical(D[Control, , drop = FALSE], Z$Death$MAR[Control, , drop = FALSE]) ||
        !identical(T[Control, , drop = FALSE], Z$End$MAR[Control, , drop = FALSE])) stop("Reference arm completions differ across scenarios.")
    Peers <- Scenarios$Scenario[Scenarios$DeathBorrowing == Scenarios$DeathBorrowing[match(Name, Scenarios$Scenario)]]
    for (Peer in Peers) if (!identical(D, Z$Death[[Peer]]) || !identical(T, Z$End[[Peer]]) ||
      !identical(Z$ExpectedDeath[[Name]], Z$ExpectedDeath[[Peer]]) || !identical(Z$ExpectedEnd[[Name]], Z$ExpectedEnd[[Peer]])) stop("Survival completions differ within a mortality rule.")
    Analyses[[Name]] <- rbind(
      GetCompletedAnalyses(C, Z$ExpectedCounts[[Name]], Target, Name, "HospitalizationCount", Meta),
      GetCompletedAnalyses(D, Z$ExpectedDeath[[Name]], Target, Name, "DeathRisk", Meta),
      GetCompletedAnalyses(T, Z$ExpectedEnd[[Name]], Target, Name, "RestrictedMeanSurvival", Meta))
    CountRows <- Analyses[[Name]]$Outcome == "HospitalizationCount"
    Active <- Target$Treatment == 1
    Analyses[[Name]]$ConditionalCountEstimate <- NA_real_
    Analyses[[Name]]$ConditionalCountEstimate[CountRows] <- colMeans(Z$ConditionalCounts[[Name]][Active, , drop = FALSE]) -
      colMeans(Z$ConditionalCounts[[Name]][!Active, , drop = FALSE])
    Burden[[Name]] <- do.call(rbind, lapply(0:1, function(Arm) {
      Keep <- Target$Treatment == Arm
      data.frame(Scenario = Name, Treatment = Arm, Participants = sum(Keep), NeedImputation = sum(Need & Keep),
        ObservedHospitalizations = sum(Target$ObservedCount[Keep]),
        MeanImputedHospitalizations = mean(colSums(Missing[Keep, , drop = FALSE])),
        MeanCompletedHospitalizations = mean(colSums(C[Keep, , drop = FALSE])),
        MeanAdditionalDeaths = mean(colSums(D[Keep, , drop = FALSE])) - sum(Target$KnownDeath[Keep]),
        MeanCompletedDeaths = mean(colSums(D[Keep, , drop = FALSE])),
        MeanAdditionalPersonYears = mean(colSums(T[Keep, , drop = FALSE] - Target$ObservedEnd[Keep])))
    }))
    Participants[[Name]] <- data.frame(Scenario = Name, Id = Target$Id, OriginalId = Target$OriginalId,
      Treatment = Target$Treatment, NeedsImputation = Need, ObservedCount = Target$ObservedCount,
      MeanCompletedCount = rowMeans(C), MeanImputedCount = rowMeans(Missing),
      ProbabilityDeathByHorizon = rowMeans(D), MeanTimeAlive = rowMeans(T))
  }
  KeepMissing <- function(X, Integer = FALSE) lapply(X, function(A) {
    B <- A[Need, , drop = FALSE]
    if (Integer && all(B <= .Machine$integer.max)) storage.mode(B) <- "integer"
    B
  })
  list(Indices = Indices, Analyses = do.call(rbind, Analyses), Burden = do.call(rbind, Burden),
    Participants = do.call(rbind, Participants), Counts = KeepMissing(Z$Counts, TRUE),
    Death = KeepMissing(Z$Death, TRUE), TimeAlive = KeepMissing(Z$End), AuditPassed = TRUE)
}

RunApplicationImputations <- function(Draws, Target, Cutpoints, Scenarios, Settings, Directory) {
  M <- CheckMiDraws(Draws, Target, Cutpoints)
  Seed <- Settings$Seed + 10000
  Batches <- split(seq_len(M), ceiling(seq_len(M)/Settings$ImputationBatchSize))
  Workers <- min(Settings$CoreBudget, length(Batches))
  Cluster <- NULL
  on.exit(if (!is.null(Cluster)) parallel::stopCluster(Cluster), add = TRUE)
  if (Workers > 1) {
    Cluster <- parallel::makePSOCKcluster(Workers)
    Functions <- c("SetSimulationSeed", "GetIntervalExposure", "DrawTerminalFollowup", "GetFutureMoments",
      "ImputeMiBatch", "GetCompletedAnalyses", "SummarizeMiBatch")
    Export <- list2env(mget(Functions, envir = environment(RunApplicationImputations)), parent = emptyenv())
    for (Name in c("Draws", "Target", "Cutpoints", "Scenarios", "Seed")) assign(Name, get(Name), Export)
    parallel::clusterExport(Cluster, c(Functions, "Draws", "Target", "Cutpoints", "Scenarios", "Seed"), envir = Export)
    parallel::clusterCall(Cluster, function() {
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      NULL
    })
  }
  Task <- function(Indices) SummarizeMiBatch(Indices, Draws, Target, Cutpoints, Scenarios, Seed)
  if (Workers > 1) environment(Task) <- .GlobalEnv
  message("Completing ", M, " datasets per scenario on ", Workers, " worker(s).")
  Results <- vector("list", length(Batches))
  for (First in seq(1, length(Batches), by = Workers)) {
    At <- First:min(length(Batches), First + Workers - 1)
    Results[At] <- if (Workers > 1) parallel::parLapply(Cluster, Batches[At], Task) else lapply(Batches[At], Task)
    if (First == 1 || First + Workers > length(Batches) || First%%max(1, 10*Workers) == 1) {
      message("  ", max(unlist(Batches[At])), "/", M, " imputations completed.")
    }
  }
  if (!identical(as.numeric(unlist(lapply(Results, `[[`, "Indices"))), as.numeric(seq_len(M))) ||
      !all(vapply(Results, `[[`, logical(1), "AuditPassed"))) stop("Imputation batches are incomplete.")
  Analyses <- do.call(rbind, lapply(Results, `[[`, "Analyses"))
  Analyses <- Analyses[order(match(Analyses$Scenario, Scenarios$Scenario), Analyses$Outcome, Analyses$Imputation), ]
  rownames(Analyses) <- NULL
  AverageTable <- function(Name, Fields) {
    X <- Results[[1]][[Name]]
    for (Field in Fields) X[[Field]] <- Reduce(`+`, lapply(Results, function(R) R[[Name]][[Field]]*length(R$Indices)))/M
    X
  }
  Burden <- AverageTable("Burden", c("MeanImputedHospitalizations", "MeanCompletedHospitalizations",
    "MeanAdditionalDeaths", "MeanCompletedDeaths", "MeanAdditionalPersonYears"))
  Burden$ImputedFractionOfMeanTotal <- with(Burden, ifelse(MeanCompletedHospitalizations > 0,
    MeanImputedHospitalizations/MeanCompletedHospitalizations, 0))
  Participants <- AverageTable("Participants", c("MeanCompletedCount", "MeanImputedCount", "ProbabilityDeathByHorizon", "MeanTimeAlive"))
  Pooled <- PoolMiAnalyses(Analyses, Settings$ConfidenceLevel)
  Last <- max(Draws$Meta$Iteration)
  Prefixes <- sort(unique(c(c(250, 500, 1000, 2000, 4000)[c(250, 500, 1000, 2000, 4000) < Last], Last)))
  Precision <- do.call(rbind, lapply(Prefixes, function(Prefix) {
    X <- PoolMiAnalyses(Analyses[Analyses$Iteration <= Prefix, ], Settings$ConfidenceLevel)
    Match <- match(paste(X$Scenario, X$Outcome), paste(Pooled$Scenario, Pooled$Outcome))
    X$IterationsPerChain <- Prefix
    X$PointMinusFinal <- X$Estimate - Pooled$Estimate[Match]
    X$SeRelativeToFinal <- X$RubinSe/Pooled$RubinSe[Match]
    X
  }))
  Shifts <- do.call(rbind, lapply(Scenarios$Scenario, function(Name) do.call(rbind, lapply(unique(Analyses$Outcome), function(Outcome) {
    A <- Analyses[Analyses$Scenario == Name & Analyses$Outcome == Outcome, ]
    B <- Analyses[Analyses$Scenario == "MAR" & Analyses$Outcome == Outcome, ]
    if (!identical(A$Imputation, B$Imputation)) stop("Paired scenario draws are misaligned.")
    Difference <- A$Estimate - B$Estimate
    data.frame(Scenario = Name, Outcome = Outcome, DifferenceFromMAR = mean(Difference),
      PairedNumericalMcse = GetChainMcse(Difference, A))
  }))))
  Combined <- function(Name) setNames(lapply(Scenarios$Scenario, function(S) do.call(cbind, lapply(Results, function(X) X[[Name]][[S]]))), Scenarios$Scenario)
  Completed <- list(Version = 1, Storage = "Only rows requiring completion; other outcomes in Target",
    Target = Target, MissingRows = which(Target$NeedsImputation), DrawIndex = Draws$Meta,
    Counts = Combined("Counts"), Death = Combined("Death"), TimeAlive = Combined("TimeAlive"), Scenarios = Scenarios, Seed = Seed)
  SaveSimulationRds(Completed, file.path(Directory, "CompletedMiData.rds"))
  Audit <- data.frame(Scenario = Scenarios$Scenario, Imputations = M, MissingParticipants = sum(Target$NeedsImputation),
    ObservedCountsUnchanged = TRUE, KnownTerminalHistoryUnchanged = TRUE, IntegerCounts = TRUE,
    CountsStopAtDeathOrHorizon = TRUE, IdenticalControlCounts = TRUE, IdenticalControlTerminal = TRUE,
    TerminalCompletionsSharedWithinRule = TRUE)
  list(Tables = list(MiAnalyses = Analyses, MiPooled = Pooled, MiBurden = Burden, MiParticipants = Participants,
    MiPrecision = Precision, MiScenarioShifts = Shifts, MiAudit = Audit), Completed = Completed)
}

ExpandCompletedOutcome <- function(Completed, Scenario = "MAR", Imputation = 1) {
  if (is.character(Completed) && length(Completed) == 1) Completed <- readRDS(Completed)
  if (!Scenario %in% names(Completed$Counts) || length(Imputation) != 1 || !is.finite(Imputation) ||
      Imputation != floor(Imputation) || Imputation < 1 || Imputation > nrow(Completed$DrawIndex)) stop("Invalid scenario or imputation number.")
  P <- Completed$Target
  P$CompletedCount <- P$ObservedCount
  P$CompletedDeath <- P$KnownDeath
  P$CompletedTimeAlive <- P$ObservedEnd
  At <- Completed$MissingRows
  P$CompletedCount[At] <- Completed$Counts[[Scenario]][, Imputation]
  P$CompletedDeath[At] <- Completed$Death[[Scenario]][, Imputation]
  P$CompletedTimeAlive[At] <- Completed$TimeAlive[[Scenario]][, Imputation]
  P
}

GetApplicationDiagnostics <- function(Fit, Draws, Prepared, Cells, Cutpoints) {
  A <- Fit$Draws
  P <- Prepared$Patients
  K <- length(Cutpoints) - 1
  N <- nrow(P)
  M <- nrow(Draws$Meta)
  Summary <- function(X) SummarizeCustomMatrix(matrix(X, nrow = dim(A)[1], ncol = dim(A)[2]))
  Parameters <- c("BetaR", "GammaR", "BetaD", "GammaD", "Alpha", "FrailtyVariance", "FrailtyPrecision",
    paste0("RecurrenceLogRate[", seq_len(K), "]"), paste0("DeathLogRate[", seq_len(K), "]"))
  ParameterSummary <- do.call(rbind, lapply(Parameters, function(Name) data.frame(Parameter = Name, t(Summary(A[, , Name])), row.names = NULL)))
  Scores <- do.call(rbind, lapply(1:5, function(J) {
    X <- A[, , paste0("Scores[", J, "]")]
    E <- posterior::mcse_mean(X)
    data.frame(Parameter = c("BetaR", "GammaR", "BetaD", "GammaD", "Alpha")[J], MeanScore = mean(X),
      Mcse = E, StandardizedScore = mean(X)/E, Review = !is.finite(E) || abs(mean(X)) > 3*E)
  }))
  Baselines <- Risk <- list()
  RateMeans <- list()
  TotalChecks <- list()
  for (Endpoint in c("Recurrence", "Death")) {
    Recurrence <- Endpoint == "Recurrence"
    LogBaseline <- if (Recurrence) Draws$RecurrenceBaseline else Draws$DeathBaseline
    Frailty <- if (Recurrence) Draws$RecurrenceFrailty else Draws$DeathFrailty
    Coeff <- Draws$Fixed[, if (Recurrence) 1:2 else 3:4, drop = FALSE]
    Eta <- Frailty + outer(P$Treatment, Coeff[, 1]) + outer(P$Age60, Coeff[, 2])
    Rates <- matrix(NA_real_, N, K)
    Totals <- numeric(M)
    for (J in seq_len(K)) {
      Base <- exp(LogBaseline[J, ])
      Baselines[[length(Baselines) + 1]] <- data.frame(Endpoint = Endpoint, Interval = J,
        Start = Cutpoints[J], Stop = Cutpoints[J + 1], t(Summary(Base)), row.names = NULL)
      Rate <- exp(sweep(Eta, 2, LogBaseline[J, ], "+"))
      if (any(!is.finite(Rate))) stop("Nonfinite posterior intensity.")
      Rates[, J] <- rowMeans(Rate)
      Exposure <- pmax(0, pmin(P$Stop, Cutpoints[J + 1]) - Cutpoints[J])
      for (Arm in 0:1) {
        Keep <- P$Treatment == Arm
        Integrated <- colSums(Rate[Keep, , drop = FALSE]*Exposure[Keep])
        Totals <- Totals + Integrated
        Observed <- sum(Cells[[if (Recurrence) "Recurrence" else "Death"]][Cells$Interval == J & Cells$Treatment == Arm])
        Risk[[length(Risk) + 1]] <- data.frame(Endpoint = Endpoint, Treatment = Arm, Interval = J,
          Start = Cutpoints[J], Stop = Cutpoints[J + 1], Observed = Observed, t(Summary(Integrated)), row.names = NULL)
      }
    }
    Stored <- as.numeric(A[, , paste0("ExpectedTotals[", if (Recurrence) 1 else 2, "]")])
    Gap <- max(abs(Totals - Stored)/(1 + abs(Stored)))
    if (!is.finite(Gap) || Gap > 1e-8) stop("Integrated intensity reconstruction differs from the sampler.")
    TotalChecks[[Endpoint]] <- data.frame(Endpoint = Endpoint, Observed = sum(if (Recurrence) P$Count else P$Death),
      t(Summary(Totals)), MaximumReconstructionDifference = Gap, row.names = NULL)
    RateMeans[[Endpoint]] <- Rates
  }
  Exposure <- GetIntervalExposure(P$Start, P$Stop, Cutpoints)
  Participant <- data.frame(P, RecurrenceExpected = rowSums(RateMeans$Recurrence*Exposure),
    DeathExpected = rowSums(RateMeans$Death*Exposure))
  Participant$RecurrenceMartingale <- Participant$Count - Participant$RecurrenceExpected
  Participant$DeathMartingale <- Participant$Death - Participant$DeathExpected
  Gaps <- do.call(rbind, lapply(seq_len(N), function(I) {
    Times <- sort(Prepared$EventData$Time[Prepared$EventData$Id == P$Id[I] & Prepared$EventData$Status == 1])
    data.frame(Id = P$Id[I], Treatment = P$Treatment[I], Age60 = P$Age60[I],
      Start = c(0, Times), Stop = c(Times, P$Stop[I]), Event = c(rep(1, length(Times)), 0))
  }))
  Gaps$Included <- Gaps$Stop > Gaps$Start
  Valid <- Gaps[Gaps$Included, ]
  GapExposure <- GetIntervalExposure(Valid$Start, Valid$Stop, Cutpoints)
  R <- rowSums(GapExposure*RateMeans$Recurrence[Valid$Id, , drop = FALSE])
  Residuals <- rbind(data.frame(Endpoint = "Recurrence", Valid, Residual = R),
    data.frame(Endpoint = "Death", Id = P$Id, Treatment = P$Treatment, Age60 = P$Age60,
      Start = 0, Stop = P$Stop, Event = P$Death, Included = TRUE, Residual = Participant$DeathExpected))
  Curves <- list()
  for (Endpoint in c("Recurrence", "Death")) for (Group in c("All", "Usual care", "Exercise")) {
    X <- Residuals[Residuals$Endpoint == Endpoint, ]
    if (Group != "All") X <- X[X$Treatment == as.numeric(Group == "Exercise"), ]
    Curves[[length(Curves) + 1]] <- data.frame(Source = "Mcmc", Endpoint = Endpoint, Group = Group,
      GetNelsonAalen(X$Residual, X$Event))
  }
  list(ParameterSummary = ParameterSummary, ScoreChecks = Scores, McmcDiagnostics = Fit$VariableDiagnostics,
    SamplingBlocks = Fit$BlockDiagnostics, ConditionalDensityAudit = Fit$DensityAudit,
    BaselineHazards = do.call(rbind, Baselines), RiskSetCalibration = do.call(rbind, Risk),
    IntegratedIntensity = do.call(rbind, TotalChecks), ResidualParticipants = Participant,
    CoxSnellResiduals = Residuals, CoxSnellCurves = do.call(rbind, Curves), GapAudit = Gaps)
}


WriteSamplerPlots <- function(Fit, Path, PlotSettings) {
  P <- CheckPlotSettings(PlotSettings)
  P$YLabelLine <- P$SamplerYLabelLine
  Names <- c("BetaR", "GammaR", "BetaD", "GammaD", "Alpha", "FrailtyVariance")
  Labels <- c("Recurrence treatment", "Recurrence age", "Mortality treatment",
    "Mortality age", "Association", "Frailty variance")
  A <- Fit$Draws
  OpenApplicationPdf(Path, P, P$SamplerFigureWidth, P$SamplerFigureHeight)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (Kind in c("Trace", "Autocorrelation")) {
    SetApplicationPanels(P, Rows = 3, Columns = 2)
    graphics::par(mai = c(0.76, 0.92, 0.12, 0.12), omi = c(P$SamplerLegendSpace, 0, 0, 0))
    for (J in seq_along(Names)) {
      X <- A[, , Names[J]]
      if (Kind == "Trace") {
        At <- unique(c(seq(1, nrow(X), by = max(1, floor(nrow(X)/2000))), nrow(X)))
        DrawApplicationPlot(P, At, X[At, 1], "Retained iteration", Labels[J], type = "n", ylim = range(X))
        for (Chain in seq_len(ncol(X))) graphics::lines(At, X[At, Chain], col = P$ArmColours[1 + (Chain - 1)%%2], lwd = 0.5)
      } else {
        Acf <- lapply(seq_len(ncol(X)), function(Chain) stats::acf(X[, Chain], lag.max = 100, plot = FALSE))
        DrawApplicationPlot(P, 0:100, as.numeric(Acf[[1]]$acf), "Lag", paste0(Labels[J], "\nAutocorrelation"),
          type = "n", ylim = c(-0.1, 1))
        graphics::abline(h = 0, lty = 3, col = P$ReferenceColour)
        for (Chain in seq_along(Acf)) graphics::lines(as.numeric(Acf[[Chain]]$lag), as.numeric(Acf[[Chain]]$acf),
          col = P$ArmColours[1 + (Chain - 1)%%2], lwd = 1)
      }
    }
    graphics::legend(graphics::grconvertX(0.5, from = "ndc", to = "user"),
      graphics::grconvertY(P$SamplerLegendSpace/(2*P$SamplerFigureHeight), from = "ndc", to = "user"),
      legend = paste("Chain", seq_len(dim(A)[2])), col = rep(P$ArmColours, length.out = dim(A)[2]),
      lty = 1, lwd = P$LineWidth, bty = "n", xpd = NA, xjust = 0.5, yjust = 0.5,
      ncol = dim(A)[2], cex = P$LegendSize/P$BaseFontSize)
  }
  invisible(Path)
}

WriteApplicationOutputs <- function(Tables, Settings, Directory, PlotSettings = HfActionPlotSettings,
                                    TableSettings = HfActionTableSettings) {
  for (Folder in c("Manuscript", "Supplement", "Plots")) dir.create(file.path(Directory, Folder), recursive = TRUE, showWarnings = FALSE)
  P <- CheckPlotSettings(PlotSettings)
  WriteMiPlots(Tables, Settings, file.path(Directory, "Plots", "HospitalizationSensitivity.pdf"),
    file.path(Directory, "Supplement", "MiDiagnostics.pdf"), P)
  WriteManuscriptModelPlots(Tables, file.path(Directory, "Plots"), P)
  WriteApplicationTables(Tables, Settings, Directory, TableSettings)
  Row <- function(X) paste0(paste(X, collapse = " & "), " \\\\")
  X <- Tables$ParameterSummary
  X <- X[X$Parameter %in% c("BetaR", "GammaR", "BetaD", "GammaD", "Alpha", "FrailtyVariance"), ]
  Labels <- c(BetaR = "Recurrence: treatment", GammaR = "Recurrence: age 60 or older", BetaD = "Mortality: treatment",
    GammaD = "Mortality: age 60 or older", Alpha = "Association parameter", FrailtyVariance = "Frailty variance")
  Rows <- vapply(seq_len(nrow(X)), function(J) Row(c(Labels[X$Parameter[J]],
    LatexNumber(unlist(X[J, c("Mean", "Sd", "Lower", "Upper")]), 3))), character(1))
  WriteLatexTable(file.path(Directory, "Manuscript", "TableModelParameters.tex"),
    "Posterior summaries for the joint recurrence and mortality model.", "tab:HfModelParameters",
    paste("l", paste(rep(LatexColumn(unlist(X[, c("Mean", "Sd", "Lower", "Upper")]), 3), 4), collapse = " ")),
    Row(c("Parameter", "{Mean}", "{SD}", "{2.5\\%}", "{97.5\\%}")), Rows,
    "Treatment coefficients compare exercise with usual care; age compares participants aged 60 or older with younger participants. Intervals are marginal posterior credible intervals. Frailty variance is summarized directly from variance draws. SD denotes posterior standard deviation.", TableSettings)
  X <- Tables$MiPooled[Tables$MiPooled$Outcome == "HospitalizationCount", ]
  X <- X[match(Settings$Scenarios$Scenario, X$Scenario), ]
  Fields <- c("Estimate", "ExpectedEstimate", "NumericalMcse", "ExpectedMcse", "PredictiveSamplingDifference", "PredictiveSamplingMcse")
  Rows <- vapply(seq_len(nrow(X)), function(J) Row(c(GetScenarioLabels(Settings$Scenarios)[J], LatexNumber(unlist(X[J, Fields]), 4))), character(1))
  WriteLatexTable(file.path(Directory, "Supplement", "TableNumericalPrecision.tex"),
    "Numerical precision of hospitalization contrasts.", "tab:HfNumericalPrecision",
    paste("l", paste(rep(LatexColumn(unlist(X[, Fields]), 4), length(Fields)), collapse = " ")),
    c(Row(c("", "{MI}", "{Integrated}", "{MCSE}", "{MCSE}", "{MI minus}", "{MCSE of}")),
      Row(c("Assumption", "{estimate}", "{estimate}", "{of MI}", "{integrated}", "{integrated}", "{difference}"))), Rows,
    "Monte Carlo standard errors (MCSEs) use the retained chains and account for serial dependence. Integration removes the random death and count completions at each posterior draw. These errors concern numerical precision and are distinct from the Rubin standard errors. Every retained draw is used.", TableSettings)
  Preamble <- c("\\usepackage{booktabs}", "\\usepackage{siunitx}", "\\usepackage{caption}", "\\usepackage{ragged2e}",
    "\\usepackage{graphicx}", "\\sisetup{group-digits=false,detect-all}", "\\captionsetup{justification=justified,singlelinecheck=false}")
  writeLines(Preamble, file.path(Directory, "LatexPreamble.tex"))
  Figure <- function(File, Caption, Label) c("\\begin{figure}[!htbp]", "\\centering",
    paste0("\\includegraphics[width=\\linewidth]{Plots/", File, "}"), paste0("\\caption{", Caption, "}"),
    paste0("\\label{", Label, "}"), "\\end{figure}")
  Horizon <- format(Settings$Horizon, trim = TRUE)
  Confidence <- format(100*Settings$ConfidenceLevel, trim = TRUE)
  FigureSpecs <- list(
    Hospitalizations = list(File = "HospitalizationSensitivity.pdf", Label = "fig:HfSensitivity", Caption = paste0(
      "Difference in mean hospitalizations through ", Horizon, " years. Points show exercise minus usual care and bars show ",
      Confidence, "\\% Rubin intervals. Counts stop at death or the horizon. All participants remain in the denominator. Table~\\ref{tab:HfScenarios} gives the borrowing assumptions.")),
    Mortality = list(File = "MortalitySensitivity.pdf", Label = "fig:HfMortality", Caption = paste0(
      "Difference in death risk through ", Horizon, " years, in percentage points. Points show exercise minus usual care and bars show ",
      Confidence, "\\% Rubin intervals. Scenarios sharing the mortality weight share their survival completions.")),
    Survival = list(File = "SurvivalSensitivity.pdf", Label = "fig:HfSurvival", Caption = paste0(
      "Difference in restricted mean survival through ", Horizon, " years, in days. Points show exercise minus usual care and bars show ", Confidence, "\\% Rubin intervals.")),
    ModelObservedFit = list(File = "ModelObservedFit.pdf", Label = "fig:SuppModelObservedFit", Caption =
      "Observed events and fitted integrated intensities over observed follow-up. Hospitalizations are above and deaths below; usual care is on the left and exercise on the right. Points and bars for fitted intensities show posterior means and 95\\% credible intervals. Death integrals are compensators over observed exposure, rather than replicated death counts."),
    ModelBaselineHazards = list(File = "ModelBaselineHazards.pdf", Label = "fig:SuppModelBaselineHazards", Caption =
      "Conditional baseline recurrence intensity and mortality hazard. Points show posterior medians and bars show 95\\% marginal credible intervals. Rates refer to usual care, age below 60, and zero log frailty. Each process has eight separate piecewise constant rates."),
    ModelCoxSnell = list(File = "ModelCoxSnell.pdf", Label = "fig:SuppModelCoxSnell", Caption =
      "Cox--Snell curves constructed from posterior mean conditional rates. The dotted identity line is the exponential reference. Curves stop when fewer than ten observations remain in the supported risk set. Positive recurrence gaps contribute event residuals and final gaps are censored. Rounded gaps of zero length are omitted from this plot only. These descriptive checks reuse fitted participant effects and observed histories.")
  )
  for (Name in names(FigureSpecs)) {
    Spec <- FigureSpecs[[Name]]
    Folder <- if (startsWith(Name, "Model")) "Supplement" else "Manuscript"
    writeLines(Figure(Spec$File, Spec$Caption, Spec$Label), file.path(Directory, Folder, paste0("Figure", Name, ".tex")))
  }
  Read <- function(Relative) readLines(file.path(Directory, Relative), warn = FALSE)
  Main <- c("TableObservedData", "TableScenarios", "TableModelParameters", "TableHospitalizations",
    "FigureHospitalizations", "TableSurvival", "FigureMortality", "FigureSurvival")
  Supp <- c("TableImputationBurden", "TableNumericalPrecision", "FigureModelObservedFit", "FigureModelBaselineHazards", "FigureModelCoxSnell")
  Preview <- c("\\documentclass[11pt]{article}", "\\usepackage[a4paper,margin=20mm]{geometry}", "\\usepackage[T1]{fontenc}",
    Preamble, "\\begin{document}", unlist(lapply(Main, function(Name) c(Read(paste0("Manuscript/", Name, ".tex")), "\\clearpage"))),
    "\\renewcommand{\\thetable}{S\\arabic{table}}", "\\setcounter{table}{0}",
    unlist(lapply(Supp, function(Name) c(Read(paste0("Supplement/", Name, ".tex")), "\\clearpage"))), "\\end{document}")
  writeLines(Preview, file.path(Directory, "ApplicationOutputs.tex"))
  Guide <- c("HF-ACTION application outputs", "",
    "ApplicationOutputs.tex and Plots contain the complete LaTeX preview.",
    "Individual table and figure fragments are in Manuscript and Supplement; LatexPreamble.tex lists the required packages.",
    "Supplement contains sampler and imputation diagnostic plots.",
    "HfActionCustomResults.rds stores results, settings, and the observed data.",
    "PosteriorDraws.rds stores all retained draws with chain and iteration indices.",
    "CompletedMiData.rds stores outcomes for participants requiring completion and the information needed to reconstruct complete datasets.",
    "WriteHfActionManuscript() creates a new export using saved results and the selected plot and table settings.")
  writeLines(Guide, file.path(Directory, "OutputGuide.txt"))
  saveRDS(list(PlotSettings = P, TableSettings = TableSettings, AnalysisSettings = Settings), file.path(Directory, "ExportSettings.rds"))
  invisible(Directory)
}

WriteHfActionManuscript <- function(Results = Application, OutputDir = NULL,
                                  PlotSettings = HfActionPlotSettings, TableSettings = HfActionTableSettings) {
  if (is.character(Results) && length(Results) == 1) Results <- readRDS(Results)
  if (is.null(Results$Tables$MiPooled) || is.null(Results$Settings)) stop("Supply results from this custom application.")
  if (is.null(OutputDir)) OutputDir <- file.path(dirname(Results$Directory), "ManuscriptExports")
  Directory <- CreateOutputDirectory(OutputDir, "Export")
  WriteApplicationOutputs(Results$Tables, Results$Settings, Directory, PlotSettings, TableSettings)
  message("Manuscript files: ", Directory)
  invisible(Directory)
}


RunHfActionCustom <- function(Data = NULL, Settings = HfActionSettings, OutputDir = NULL,
                             Scenarios = GetMiScenarios(), PlotSettings = HfActionPlotSettings,
                             TableSettings = HfActionTableSettings) {
  Started <- proc.time()[3]
  ValidateApplicationSettings(Settings)
  PlotSettings <- CheckPlotSettings(PlotSettings)
  Scenarios <- CheckMiScenarios(Scenarios)
  Software <- EnsureCustomPackages(Settings)
  if (!Settings$InstallMissingPackages && is.null(Data) && !requireNamespace("rmt", quietly = TRUE)) stop("Install rmt to load the public dataset.")
  if (is.null(OutputDir)) {
    if (!nzchar(HfActionProgramPath)) stop("Run this script with source() or Rscript.")
    OutputDir <- file.path(dirname(HfActionProgramPath), Settings$OutputFolder)
  }
  Prepared <- ReadData(Data)
  P <- Prepared$Patients
  Target <- PrepareMiTarget(Prepared, Settings$Horizon)
  Cutpoints <- Settings$Cutpoints
  if (is.null(Cutpoints)) Cutpoints <- seq(0, max(P$Stop), length.out = Settings$Intervals + 1)
  if (!is.numeric(Cutpoints) || length(Cutpoints) < 3 || any(!is.finite(Cutpoints)) || Cutpoints[1] != 0 ||
      any(diff(Cutpoints) <= 0) || tail(Cutpoints, 1) < max(P$Stop)) stop("Invalid baseline cutpoints.")
  Cells <- MakeRiskCells(Prepared, Cutpoints)
  for (I in seq_len(nrow(P))) {
    C <- Cells[Cells$Id == P$Id[I], ]
    if (sum(C$Recurrence) != P$Count[I] || sum(C$Death) != P$Death[I] || abs(sum(C$Exposure) - P$Stop[I]) > 1e-10) stop("Participant likelihood audit failed.")
  }
  Cores <- parallel::detectCores(logical = TRUE)
  if (length(Cores) != 1 || !is.finite(Cores) || Cores < 1) Cores <- 1
  Settings$LogicalCores <- Cores
  Settings$CoreBudget <- max(1, floor(Settings$CoreFraction*Cores))
  Settings$Cutpoints <- Cutpoints
  Settings$Intervals <- length(Cutpoints) - 1
  Settings$Scenarios <- Scenarios
  Settings$Version <- 1
  Settings$Software <- Software
  Settings$TimeUnit <- "years"
  Directory <- CreateOutputDirectory(OutputDir, "Run")
  on.exit(message("Run folder: ", Directory), add = TRUE)
  writeLines(Directory, file.path(OutputDir, "LatestApplicationDirectory.txt"))
  if (nzchar(HfActionProgramPath) && file.exists(HfActionProgramPath)) file.copy(HfActionProgramPath,
    file.path(Directory, "ProgramUsed.R"), overwrite = FALSE)
  SaveSimulationRds(list(Settings = Settings, Prepared = Prepared, Target = Target, Cells = Cells), file.path(Directory, "PreparedInput.rds"))
  message("HF-ACTION custom application: ", nrow(P), " participants, ", sum(P$Count), " hospitalizations, ",
    sum(P$Death), " deaths; ", Settings$Intervals, " baseline intervals.")
  Fit <- FitApplication(Prepared, Cells, Cutpoints, Settings, Directory)
  Settings$Imputations <- prod(dim(Fit$Draws)[1:2])
  Settings$RetainedPerChain <- Fit$DrawsPerChain
  Settings$ActualWarmup <- Fit$WarmupIterations
  Settings$KernelHash <- Fit$ConditionalKernelHash
  Draws <- ExtractApplicationDraws(Fit, Target, Cutpoints)
  SaveSimulationRds(list(Settings = Settings, Draws = Fit$Draws, DrawIndex = Draws$Meta,
    SamplingBlocks = Fit$BlockDiagnostics, ChainSummary = Fit$Runs), file.path(Directory, "PosteriorDraws.rds"))
  message("Calculating posterior summaries and residual diagnostics...")
  Tables <- GetApplicationDiagnostics(Fit, Draws, Prepared, Cells, Cutpoints)
  Tables$ArmSummary <- do.call(rbind, lapply(0:1, function(Arm) {
    X <- P[P$Treatment == Arm, ]
    data.frame(Treatment = Arm, Participants = nrow(X), Hospitalizations = sum(X$Count), Deaths = sum(X$Death),
      Censored = sum(X$Death == 0), PersonYears = sum(X$Stop))
  }))
  Tables$HorizonSummary <- do.call(rbind, lapply(0:1, function(Arm) {
    X <- Target[Target$Treatment == Arm, ]
    data.frame(Treatment = Arm, Horizon = Settings$Horizon, Participants = nrow(X), NeedImputation = sum(X$NeedsImputation),
      PercentNeedImputation = 100*mean(X$NeedsImputation), ObservedDeathsByHorizon = sum(X$KnownDeath),
      ObservedThroughHorizon = sum(X$ObservedEnd == Settings$Horizon & X$KnownDeath == 0),
      ObservedHospitalizationsByHorizon = sum(X$ObservedCount), ObservedPersonYears = sum(X$ObservedEnd))
  }))
  Tables$RiskCells <- Cells
  Tables$BoundaryEvents <- Prepared$Boundary
  Tables$MiScenarios <- Scenarios
  dir.create(file.path(Directory, "Supplement"), showWarnings = FALSE)
  WriteSamplerPlots(Fit, file.path(Directory, "Supplement", "SamplerDiagnostics.pdf"), PlotSettings)
  FitSummary <- Fit[setdiff(names(Fit), c("Draws", "Data", "VariableDiagnostics", "DensityAudit"))]
  rm(Fit)
  invisible(gc())
  Mi <- RunApplicationImputations(Draws, Target, Cutpoints, Scenarios, Settings, Directory)
  Tables <- c(Tables, Mi$Tables)
  rm(Mi, Draws)
  invisible(gc())
  Settings$TotalSeconds <- unname(proc.time()[3] - Started)
  Result <- list(Version = 1, Directory = Directory, Settings = Settings, Tables = Tables, FitSummary = FitSummary,
    RawData = Prepared$RawData, Patients = Prepared$Patients, Target = Target)
  SaveSimulationRds(Result, file.path(Directory, "HfActionCustomResults.rds"))
  for (Name in names(Tables)) if (is.data.frame(Tables[[Name]])) utils::write.csv(Tables[[Name]],
    file.path(Directory, paste0(Name, ".csv")), row.names = FALSE, na = "")
  WriteApplicationOutputs(Tables, Settings, Directory, PlotSettings, TableSettings)
  Pooled <- Tables$MiPooled
  Diagnostics <- Tables$McmcDiagnostics
  Text <- c("HF-ACTION application", "",
    paste("Participants:", nrow(P), "Observed hospitalizations:", sum(P$Count), "Observed deaths:", sum(P$Death)),
    paste("Horizon (years):", Settings$Horizon, "Participants requiring completion:", sum(Target$NeedsImputation)),
    paste("Separate baseline intervals per process:", Settings$Intervals),
    paste("Fixed effect SD:", Settings$FixedEffectSd, "Alpha SD:", Settings$AlphaSd,
      "Baseline log mean/SD:", Settings$BaselineLogMean, Settings$BaselineLogSd),
    paste("Frailty precision Gamma shape/rate:", Settings$FrailtyShape, Settings$FrailtyRate),
    paste("Chains:", Settings$Chains, "Warmup per chain:", FitSummary$WarmupIterations,
      "Retained per chain:", Settings$RetainedPerChain, "Stop:", FitSummary$StopReason),
    paste("Maximum Rhat:", signif(max(Diagnostics$Rhat), 5), "Minimum bulk/tail ESS:",
      signif(min(Diagnostics$BulkEss), 5), signif(min(Diagnostics$TailEss), 5)),
    paste("MCMC flags:", sum(Diagnostics$Review), "Score flags:", sum(Tables$ScoreChecks$Review),
      "Sampler warnings:", length(FitSummary$Warnings)),
    if (length(FitSummary$Warnings)) paste("Sampler warning:", FitSummary$Warnings),
    paste("Sampling wall seconds:", round(FitSummary$SamplingSeconds, 2), "Analysis seconds before final exports:", round(Settings$TotalSeconds, 2)),
    paste("Imputations per scenario:", Settings$Imputations),
    "Monte Carlo standard errors account for dependence within chains.", "",
    "Posterior parameter summaries:", capture.output(print(Tables$ParameterSummary, row.names = FALSE)), "",
    "Pooled results (death risk in proportions, restricted mean survival in years):",
    capture.output(print(Pooled[, c("Scenario", "Outcome", "Estimate", "RubinSe", "Lower", "Upper", "NumericalMcse", "NumericalVarianceFraction")], row.names = FALSE)),
    "", "Paired scenario changes:", capture.output(print(Tables$MiScenarioShifts, row.names = FALSE)), "",
    "Automatic sampling blocks:", capture.output(print(Tables$SamplingBlocks, row.names = FALSE)), "",
    "Posterior score checks:", capture.output(print(Tables$ScoreChecks, row.names = FALSE)))
  writeLines(Text, file.path(Directory, "Summary.txt"))
  message("Results: ", Directory)
  if (any(Diagnostics$Review) || any(Tables$ScoreChecks$Review) || length(FitSummary$Warnings)) message("Some diagnostics require review; finite results have been retained and flagged.")
  invisible(Result)
}

if (isTRUE(getOption("HfActionAutoRun", TRUE))) {
  Application <- RunHfActionCustom()
}