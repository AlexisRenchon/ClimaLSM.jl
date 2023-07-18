# Microbial respiration
This section describes multiple models of soil organic decomposition
by microbes, implemented in ClimaLSM. 

## Dual Arrhenius Michaelis-Menten

```@raw html
<iframe src="https://clima.westus3.cloudapp.azure.com/jsserve/Rh"
   style="height:2400px;width:100%;">
</iframe>
```

The Dual Arrhenius and Michaelis-Menten (DAMM) kinetics model in ClimaLSM.jl follows Davidson et al. 2012. DAMM models heterotrophic respiration ($Rh$) as a function of soil temperature ($T_s$) and soil moisture ($\theta$). \newline

The rate of respiration, $Rh$, is expressed as:
```math
\begin{equation}
    Rh = V_\text{maxs_x}MM_{s_x}MM_{O_2}
\end{equation}
```

where $V_\text{max}{s_x}$ is the maximum potential rate of respiration, $MM_{s_x}$ represents the availability of substrate, and $MM_{O_2}$ is used as the oxygen limitation factor. $MM_{s_x}$ and $MM_{O_2}$ are between 0 (limiting) and 1 (non limiting). $V_\text{max}{s_x}$, $MM_{s_x}$, and $MM_{O_2}$ are expressed as:

```math
\begin{equation}
    V_\text{max}{s_x} = \alpha_{s_x} \exp(\frac{-Ea_{s_x}}{RT_s})
\end{equation}
```
```math
\begin{equation}
    MM_\text{sx} = \frac{[s_x]}{kM_{s_x}+[s_x]}
\end{equation}
```
```math
\begin{equation}
    MM_\text{O_2} = \frac{[O_2]}{kM_{O_2}+[O_2]}
\end{equation}
```
where $\alpha_{s_x}$ is the pre-exponential factor, $Ea_{s_x}$ is the activation energy of the reaction, $R$ is the gas constant, and $T_s$ is soil temperature. $[s_x]$ is the concentration of all soluble substrate, and $[O_2]$ is the oxygen concentration. $kM_{s_x}$ and $kM_{O_2}$ are the Michaelis constant for soil and oxygen, respectively. \newline

The concentration of soluble carbon substrates is affected by soil water content, and specifically by diffusion of substrates through soil water films. Using these underlying principles, $[s_x]$ is calculated as:
```math
\begin{equation}
    [s_x] = p_{s_x}\times[C_{som}]\times D_{liq}\times\theta^3
\end{equation}
```
where $[C_{som}]$ is the total amount of soil organic carbon, and $p_{s_x}$ is the fraction of $[C_{som}]$ that is soluble. $D_{liq}$ is the diffusion coefficient of the soluble carbon. $\theta$ is soil moisture.

The concentration of $O_2$ depends on the diffusion of gases within the soil, which is calculated as below:
```math
\begin{equation}
    [O_2] = D_{Oa}\times O_{2a} \times porosity_{air}^{4/3}
\end{equation}
```
where $D_{Oa}$ is the diffusion coefficient for $O_2$ in air, $O_{2a}$ is the volume fraction of $O_2$ in air, and $porosity_{air}$ is the air-filled porosity.

The air-filled porosity is calculated by subtracting the soil moisture from the total porosity ($\nu$):
```math
\begin{equation}
    porosity_{air} = \nu - \theta
\end{equation}
```

To sum up, the model has the following parameters:
\baselineskip=10bp
\smallskipamount=\baselineskip
\medskipamount=2\baselineskip
\setbox\strutbox=\hbox{%
  \vrule height .7\baselineskip depth .3\baselineskip width 0pt}

\newcount\rowcount

\def\headersfor#1{
  \noalign{\global\rowcount=0 \medbreak}
  \bf #1& Symbol& Unit& Range\crcr
  \noalign{\nobreak\smallskip}}

\def\cr{\crcr\noalign{\maybeskip}}

\def\maybeskip{\ifnum\rowcount=2 \global\rowcount=0 \smallbreak
  \else \global\advance\rowcount by 1 \fi}

\halign{#\hfil\strut&& \quad\hfil#\crcr
  \headersfor{Output}
  Heterotrophic respiration&         Rh& $\mu$mol $m^{-2}$ $s^{-2}$& 0--25 \cr
  \headersfor{Drivers}
  Soil temperature& $T_s$& $\degree C$& -20--50 \cr
  Soil moisture&            $\theta$& $m^3$ $m^{-3}$& 0.0--1.0 \cr
  \headersfor{Parameters}
  Soil porosity&   $\nu$& $m^3$ $m^{-3}$& 0.0--1.0\cr
  Pre-exponential factor&  $\alpha_{s_x}$& kg C $m^{-3}$ $s^{-1}$& 100e3--300e3\cr
  Activation energy& $Ea_{s_x}$& J$mol^{-1}$& 50e3--70e3\cr
  Michaelis constant for soil& $kM_{s_x}$& kg C $m^{-3}$& 1e-10--0.1\cr
  Michaelis constant for $O_2$& $kM_{O_2}$& $m^3$$$m^{-3}$& 1e-10--0.1\cr
  Volumetric fraction of $O_2$ in the soil air content&      $O_{2_a}$& -& 0.005--0.5\cr
  Fraction of soil carbon that is considered soluble&   $p_{s_x}$& -& 0.005--0.5\cr
  Soil organic C&   $C_{som}$& kg C $m^{-3}$& 1.0--10.0\cr
}

\def\headersfor#1{
  \noalign{\global\rowcount=0 \medbreak}
  \bf #1& Symbol& Unit& Value\crcr
  \noalign{\nobreak\smallskip}}

\def\cr{\crcr\noalign{\maybeskip}}

\def\maybeskip{\ifnum\rowcount=2 \global\rowcount=0 \smallbreak
  \else \global\advance\rowcount by 1 \fi}

\halign{#\hfil\strut&& \quad\hfil#\crcr  
  \headersfor{Constants}
  Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 Pa)& $O_{a100}$& -& 0.1816\cr
  Diffusivity of soil C substrate in liquid&    $D_{liq}$& -& 3.17\cr
  Diffusion coefficient of oxygen in air&          $D_{Oa}$& -& 1.67\cr
}