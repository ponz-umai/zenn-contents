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



これは水素様原子中の電子のシュレーディンガー方程式

$$
\left(-\frac{\hbar^2}{2m} \nabla^2  -\frac{Ze^2}{4\pi\epsilon _0r}  \right) \varphi(\boldsymbol{r}) =\epsilon \varphi (\boldsymbol{r} )
$$

において第2項に$V(r)$が付け加わった形になっています。

したがって水素原子の場合で行ったように極座標表示に変換して変数分離をすることで、3つの微分方程式

$$
-\frac{\hbar^2}{2m} \left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        -\frac{\lambda }{r^2} R
        \right)
        +
        \left(
        -\frac{Ze^2}{4\pi\epsilon _0r}
        + V(r)
        \right)R = \epsilon R
       ,\\
\frac{1}{\sin\theta } \frac{d}{d\theta } \left( \sin\theta \frac{d\Theta }{d\theta }   \right) + \left( \lambda - \frac{m^2}{\sin^2\theta }  \right) \Theta =0,\\
\frac{d^2\Phi }{d\phi ^2} + m^2\Phi =0.
$$

と、クーロンポテンシャル部分に新たに平均場ポテンシャル$V(r)$が付け加わった形に整理できます。

ところで今回、一体の固有関数は球対称と仮定したのでした。そこで角度部分は定数、具体的には
$$
Y_0^0 = \frac{1}{\sqrt{4\pi}}
$$

となり、
試行関数$\varphi(|\boldsymbol{r}|)$を求めるための微分方程式は、動径部分のみの1変数微分方程式で$\lambda = 0$とした
