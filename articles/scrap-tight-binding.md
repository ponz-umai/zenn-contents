---
title: "没原稿置き場"
emoji: "👏"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: []
published: false
---
それには束縛条件の下での最適化法「ラグランジュの未定乗数法」を用いて、束縛条件は

$$
\begin{align*}
&\int\Phi^*(\tau_1,\tau_2)\Phi(\tau_1,\tau_2) d\tau_1d\tau_2\\
&= \int \varphi^*(1)\varphi^*(2)\varphi(1)\varphi(2)d\boldsymbol{r}_1d\boldsymbol{r}_2\\
&=\int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}\int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}\\
&=1
\end{align*}
$$

となり結局$\int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}=1$を満たせばよいので、未定乗数$\varepsilon$を用いて、次の式を最小にする関数$\varphi$を見つければよいことになります。

$$
\int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1d\boldsymbol{r}_2 - \varepsilon
\left[
      \int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}-1  
\right]


$$
------------------------
$$
\begin{align*}
&=
 2\int
    \delta\varphi^*(\boldsymbol{r})
    \hat{H}
    \varphi(\boldsymbol{r})
    d\boldsymbol{r}
    +
    2\int
        \varphi^*(\boldsymbol{r})
    \hat{H}
    \delta\varphi(\boldsymbol{r})
    d\boldsymbol{r}
\end{align*}
$$