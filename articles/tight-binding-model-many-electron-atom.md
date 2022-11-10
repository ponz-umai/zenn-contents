---
title: "多電子原子内の電子の波動関数と平均場ポテンシャル"
emoji: "🌌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: false
---
# はじめに
[前回記事](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)でスピン状態まで考慮した多電子の波動関数について理解しました。

前回記事では、まず一体のポテンシャルのみを取り入れて多電子系の波動関数の性質を理解しました。

その内容を踏まえて、いよいよ今回は多電子原子中の電子について取り扱います。

# He原子

まずはもっとも簡単な多電子原子、He原子のシュレディンガー方程式を題材に、多電子原子の波動関数を理解していきます。例によって雑なイラストで表すと、

![](images/tb/he.png)

こんな感じです。

[水素様原子の場合](https://zenn.dev/ponzumai/articles/tight-binding-model-hydrogen-atom)は、重心運動と相対運動に分離することにより、原子核位置を原点とした電子の相対座標のシュレディンガー方程式を導出しました。
しかし今回は3体の座標を考える必要があります。個の場合にきれいに運動を分離できるかどうか私は知りません。

3体以上が絡む運動を考える場合、思い切って原子核が静止していることにして問題を解いていくことにします。これを**ボルンーオッペンハイマー近似（ Born–Oppenheimer approximation）**と言います。大層なことを書きましたが、とにかく原子核の運動は考えないことにして、原子核位置を原点とし、2電子の運動のみを考えていくことにします。

そこで2電子の従うシュレディンガー方程式は、前章と同様にスピン座標も含めた座標を$\tau = (\boldsymbol{r},\sigma)$として以下のようになります。

$$
\mathcal{H} \Phi_{12}(\tau_1,\tau_2) = E\Phi_{12}(\tau_1,\tau_2),\\
\mathcal{H} = \frac{-\hbar^2}{2m}\nabla_1{}^2 + \frac{-\hbar^2}{2m}\nabla_2{}^2 
+ \frac{-Ze^2}{4\pi\varepsilon_0r_1}
+ \frac{-Ze^2}{4\pi\varepsilon_0r_2}
+ \frac{e^2}{4\pi\varepsilon_0r_{12}}
$$

ここで$\mathcal{H}$の第1項・2項は運動エネルギー、第3項・4項は原子核から受ける引力クーロン相互作用ポテンシャル、第5項は電子間の斥力クーロン相互作用ポテンシャルで、今ヘリウム原子を考えているので第3項・4項の原子核の電荷は$Z=2$です。

2電子のみの運動を考えると言っても、ついに2体相互作用が出てきてしまいました。一般の多電子原子も含め、この2体相互作用の扱い方が本章のメインテーマです。

## 第ゼロ近似～電子間相互作用を無視～

とはいえいきなり飛ばしていくのはよくないので、前章の復習も兼ねてまずは2体相互作用項を無視するところから始めてみます。

これだけでも結構面白いので、早速やっていきましょう。

第5項のみを無視するとシュレディンガー方程式は

$$
\left( \hat{H}_1 + \hat{H}_2 \right)\Phi_{12}(\tau_1,\tau_2) = E\Phi_{12}(\tau_1,\tau_2),\\
\hat{H}_i = \frac{-\hbar^2}{2m}\nabla_i{}^2 
+ \frac{-Ze^2}{4\pi\varepsilon_0r_i}
$$

と変数分離型で書くことができ、一体のハミルトニアンは水素様原子のハミルトニアンと同じものになります。従って一体固有関数は、水素様原子中の電子の固有関数

$$
\varphi_{nlm}(\boldsymbol{r}) = R_{nl}(r)Y_l^m(\theta,\phi)
$$

を軌道部分とし、スピン関数$\alpha(\sigma), \beta(\sigma)$をかけたものとなります。

一体の固有エネルギーは主量子数$n$のみに依存し

$$
\epsilon _n = -\frac{Z^2me^4}{(4\pi\epsilon )^22\hbar^2} \frac{1}{n^2} , \>\>\>\> n=1,2,\cdots
$$

となります。

特に基底状態の多電子波動関数は、一体の固有エネルギーが小さい順からスレーター行列式に"詰めて"行けばよかったので、$n=1$の状態$\varphi_{100}(\boldsymbol{r})$を選んで、

$$
\begin{align*}
\Phi_{12}^{GS}(\tau_1,\tau_2)  &= 
\left|
    \varphi_{100} \overline{\varphi}_{100}
\right|\\
&=\frac{1}{\sqrt{2}}\left\{
    \varphi_{100}(\boldsymbol{r}_1)\alpha(\sigma_1) 
    \varphi_{100}(\boldsymbol{r}_2)\beta(\sigma_2) 
    -
    \varphi_{100}(\boldsymbol{r}_1)\beta(\sigma_1) 
    \varphi_{100}(\boldsymbol{r}_2)\alpha(\sigma_2) 
\right\}
\end{align*}
$$

となります。

## Self-consistent法

上記電子間相互作用を無視した解を求めましたが、実際には電子はそれぞれもう一方の電子からの斥力ポテンシャルを感じているはずです。
特に基底状態の場合に対して、その効果を取り入れた近似的な波動関数の形を求める方法が**self-consistent法**（自己無撞着法）です。これは後に出てくる、一般の多電子状態に適用するHartree-Fock近似の最も簡単な場合です。

<!--参考：http://cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/atom.html-->

近似的な波動関数を求めるために、量子力学の「変分原理」の考え方（変分法）を用います。
すなわち、試行関数を用意して、ハミルトニアン演算子の期待値、つまりエネルギー期待値をを最小にするという条件を満たす関数形を求めます。

具体的には、まずは相互作用を無視した場合と同様に、基底状態の波動関数の解の形を何らかの軌道関数$\varphi(\boldsymbol{r})$にスピン関数$\alpha(\sigma), \beta(\sigma)$をそれぞれかけ合わせた関数からなるスレーター行列式と仮定します。この時軌道関数は規格化されているものと仮定します。

すなわち

$$
\begin{align*}
\Phi(\tau_1,\tau_2)  &=
\left|
    \varphi \overline{\varphi}
\right|\\
&=\frac{1}{\sqrt{2}}\left\{
    \varphi(\boldsymbol{r}_1)\alpha(\sigma_1) 
    \varphi(\boldsymbol{r}_2)\beta(\sigma_2) 
    -
    \varphi(\boldsymbol{r}_1)\beta(\sigma_1) 
    \varphi(\boldsymbol{r}_2)\alpha(\sigma_2) 
\right\}
\end{align*}
$$

ハミルトニアン演算子

$$
\mathcal{H} = \frac{-\hbar^2}{2m}\nabla_1{}^2 + \frac{-\hbar^2}{2m}\nabla_2{}^2 
+ \frac{-Ze^2}{4\pi\varepsilon_0r_1}
+ \frac{-Ze^2}{4\pi\varepsilon_0r_2}
+ \frac{e^2}{4\pi\varepsilon_0r_{12}}
$$

に対するエネルギー期待値は、途中、$\varphi(\boldsymbol{r_1}), \alpha(\sigma_1)$を$\varphi(1)\alpha(1)$等と略記して、

$$
\begin{align*}
\langle E \rangle &= \int\Phi^*(\tau_1,\tau_2) \mathcal{H}\Phi(\tau_1,\tau_2) d\tau_1d\tau_2\\
&=
    \frac{1}{2}
    \int \left[\varphi(1)\alpha(1)\varphi(2)\beta(2) -\varphi(1)\beta(1)\varphi(2)\alpha(2)\right]^*\mathcal{H}\\
    &\>\>\>\>\>\>\>\>\left[\varphi(1)\alpha(1)\varphi(2)\beta(2) -\varphi(1)\beta(1)\varphi(2)\alpha(2)\right]d\boldsymbol{r}_1d\boldsymbol{r}_2d\sigma_1\sigma_2\\
&=
    \frac{1}{2}\left[ 
    \int\alpha^*(1)\beta^*(2)\alpha(1)\beta(2)d\sigma_1d\sigma_2
    \int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1\boldsymbol{r}_2\right.\\
    &\>\>\>\>\>\>-
    \int\alpha^*(1)\beta^*(2)\beta(1)\alpha(2)d\sigma_1d\sigma_2
    \int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1\boldsymbol{r}_2\\
    &\>\>\>\>\>\>-
    \int\alpha^*(1)\beta^*(2)\beta(1)\alpha(2)d\sigma_1d\sigma_2
    \int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1\boldsymbol{r}_2\\
    &\>\>\>\>\>\>+ \left.
    \int\beta^*(1)\alpha^*(2)\beta(1)\alpha(2)d\sigma_1d\sigma_2
    \int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1\boldsymbol{r}_2
   \right]\\
&=\int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1d\boldsymbol{r}_2
\end{align*}
$$

となります。ここで2式目の第2項、第3項はスピン関数の直交性よりゼロになり、また二つの固有関数の軌道関数を同一の関数と仮定したことから、最後の式のように一つの積分だけが残ります。

ハミルトニアン演算子の中身を具体的に書いてもう少し計算を進めると、

$$
\begin{align*}
&\int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1d\boldsymbol{r}_2\\
&=
    2\int\varphi^*(\boldsymbol{r})\left( \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-Ze^2}{4\pi\varepsilon_0r} \right)
    \varphi(\boldsymbol{r})d\boldsymbol{r}\\
&\>\>\>\>+
    \int\int \varphi^*(\boldsymbol{r}_1)\varphi^*(\boldsymbol{r}_2)\frac{e^2}{4\pi\varepsilon_0\left|\boldsymbol{r}_{1}-\boldsymbol{r}_2\right|}
    \varphi(\boldsymbol{r}_1)\varphi(\boldsymbol{r}_2)d\boldsymbol{r}_1d\boldsymbol{r}_2
\end{align*}
$$

ここに、束縛条件として軌道関数の規格化条件を課した上で、上記エネルギー期待値を最小にする関数$\varphi(\boldsymbol{r})$を求めればよいことになります。

を最小にする

# 一般の多電子原子
## Hartree近似

## Hartree-Fock近似

## 多電子原子内の電子の波動関数と平均場ポテンシャル

