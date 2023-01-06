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

面心立方格子の結晶内のポテンシャル

面心立方格子は以下の基本併進ベクトルで表され、

$$
\boldsymbol{a}_1 = \frac{a}{2}\boldsymbol{e}_x + \frac{a}{2}\boldsymbol{e}_y, \\ 
{}\\
\boldsymbol{a}_2 =  \frac{a}{2}\boldsymbol{e}_y + \frac{a}{2}\boldsymbol{e}_z,\\
{}\\
\boldsymbol{a}_3 =  \frac{a}{2}\boldsymbol{e}_x + \frac{a}{2}\boldsymbol{e}_z.
$$

結晶内のポテンシャルは格子ベクトルを

$$
\boldsymbol{n} = n_1\boldsymbol{a}_1 
+
n_2\boldsymbol{a}_2
+
n_3\boldsymbol{a}_3
$$

として$V(\boldsymbol{r}+\boldsymbol{n}) = V(\boldsymbol{r})$の周期を持ちます。この時逆格子ベクトルを

$$
\boldsymbol{K} = m_1\boldsymbol{b}_1 + m_2\boldsymbol{b}_2 + m_3\boldsymbol{b}_3,\\

\boldsymbol{b}_1 
=
\frac{\boldsymbol{a}_2\times\boldsymbol{a}_3}{\boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)}
=
\frac{2\pi}{a}(\boldsymbol{e}_x +\boldsymbol{e}_y - \boldsymbol{e}_z ) 
$$

と置くと

結晶内のポテンシャルは


$$
V(\boldsymbol{r}) = \sum_{\boldsymbol{K}}V_{\boldsymbol{K}}e^{i\boldsymbol{K}\cdot\boldsymbol{r}}\\
V_{\boldsymbol{K}}= \frac{1}{v_c}\int_{v_c} V(\boldsymbol{r})e^{-i\boldsymbol{K}\cdot\boldsymbol{r}}d\boldsymbol{r}
$$


と書けます。


$H_q$を対角化する行列$T$：

$$
T^{-1}H_qT = 
\begin{bmatrix}
\ddots &   &  &  &   \\
 & \>\varepsilon_{n,q} & \\
  &&\ddots  & & &\\
 & & & \>\>\varepsilon_{n',q} \\
 & &  &  & \ddots &  &  \\
   & & &&& \varepsilon_{n'',q} & \\
&& & & &&\ddots &  \\
\end{bmatrix}
$$


を考えると（ラベルのセンスが悪い気がしますが））、「固有ベクトル」

$$
T
\begin{bmatrix}
\vdots \\
c_{q - K_l} \\
\vdots \\
c_q \\
\vdots \\
c_{q + K_l} \\
\vdots
\end{bmatrix}

=
\begin{bmatrix}
\vdots \\
\sum_{m}T_{nm}c_{q - K_m} \\
\vdots \\
\sum_{m}T_{n'm}c_{q - K_m} \\
\vdots \\
\sum_{m}T_{n''m}c_{q - K_m} \\
\vdots
\end{bmatrix}
$$

と、（無限個の）$H_{q}$に対する$n$でラベルされた固有値

$$
\varepsilon_{n,q}
$$

を得ることができます。

（行列表示の都合で$K$を$K_m$、全ての$K$に関する和を$\sum_m$と書いていたところを、2行目で「全ての$K$に関する和$\sum_K$と書き換えました）


波数空間における周期性とワニエ関数
これは次の章に持っていくか



あるいは、もう少し物理的に考えてみると、
いわゆる「空格子近似」

逆格子ベクトルだけ異なる固有エネルギーは、波数の刻みに対して離散的だとみなせるので

ーーー

元々我々は多電子の固有状態を求めたいのでした。
そこでスレーター行列式を考えなおすと、固有状態
を「詰めて」行けば多電子の固有状態を得られるわけです。

今考えている周期的固有条件に対応して、N個の格子が存在していると考えられます。
格子一つ当り1個の電子がいると考えると、ブリルアンゾーン内にはN個の固有状態が並んでおり、その一つ一つの点に対して一つの1電子固有状態が対応します。

このような考えのもとエネルギー固有値が低い順からN個、1電子固有状態を詰めていったら、最高エネルギー固有値の値を考えることができます。そのようなエネルギーをフェルミエネルギーと呼び、$\varepsilon_F$と書きます。
この値はエネルギーバンド図では等エネルギー線として横軸と平行な直線で表されます。

バンド理論超概略

スレーター行列式と一粒子固有状態
フェルミエネルギー
エネルギーギャップ
　NFEモデルー摂動展開




Wannier関数参考文献：
https://cfm.ehu.es/ivo/publications/marzari-psik03.pdf
https://www.diva-portal.org/smash/get/diva2:631027/FULLTEXT01.pdf
http://www.ms.osakafu-u.ac.jp/~taguchi/densi/5th.pdf