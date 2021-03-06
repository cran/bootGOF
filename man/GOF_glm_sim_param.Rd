% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GOF_model_simulator.R
\name{GOF_glm_sim_param}
\alias{GOF_glm_sim_param}
\title{Implements the "interface" GOF_model_simulator for
  for generalized linear models}
\description{
after the GLM was fitted the distribution of the
  of the dependent variable is fully specified and used here
  to generate new depenedent variables that follow \code{model}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-resample_y}{\code{GOF_glm_sim_param$resample_y()}}
\item \href{#method-clone}{\code{GOF_glm_sim_param$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-resample_y"></a>}}
\if{latex}{\out{\hypertarget{method-resample_y}{}}}
\subsection{Method \code{resample_y()}}{
see \link{GOF_model_simulator}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GOF_glm_sim_param$resample_y(model)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{model}}{see \link{GOF_model_simulator}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
see \link{GOF_model_simulator}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GOF_glm_sim_param$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
