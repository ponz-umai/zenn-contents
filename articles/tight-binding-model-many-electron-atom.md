---
title: "多電子原子内の電子の波動関数と平均場ポテンシャル"
emoji: "🌌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: false
---
# はじめに
[前回記事](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)でまず一体のポテンシャルのみを取り入れて、スピン状態まで考慮した多電子の波動関数について理解しました。

とはいえ、固体内の状態は（原子核が静止していると近似したとしても）10の二十数乗個の電子が相互作用しあいながら存在している状態で、到底そのようなシュレーディンガー方程式は解けません。

そこでtight-bindingモデルは、固体内の電子の状態を

- 価電子だけに注目し、内殻電子との相互作用は「価電子が感じるポテンシャル」と近似して扱う
- 価電子間のクーロン相互作用はひとまず無視する

として、一体のハミルトニアンの総和の形に近似する考え方です。

一つの価電子のみに着目して図にするとこんな感じ（例によってぐちゃぐちゃですが）

![](/images/tb/solid.png)

ハミルトニアンにすると、位置$\boldsymbol{R}$の原子核＋内殻電子から、$i$番目の電子への有効ポテンシャルを$V_i(\boldsymbol{r}_i-\boldsymbol{R})$として、

$$
\mathcal{H}\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) = E\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots),\\
\mathcal{H} = \sum_i\left(  -\frac{\hbar^2}{2m}\nabla_i^2\right) + \sum_i\sum_{\boldsymbol{R}}V_i(\boldsymbol{r}_i-\boldsymbol{R}) 
$$

となります。このように考えれば、あとは変数分離をして一電子シュレーディンガー方程式

$$
\left\{
    -\frac{\hbar^2}{2m}\nabla_i^2
    +
    \sum_{\boldsymbol{R}}V_i(\boldsymbol{r}_i-\boldsymbol{R})
\right\}
\varphi(\boldsymbol{r}) = \varepsilon\varphi(\boldsymbol{r})
$$

を解き、そこから求めた一体の固有関数をエネルギーが低い順からスレーター行列式に詰めて行けば基底状態の波動関数を得ることができます。

というわけで上記のような形をした波動関数を解いていきたいところなのですが、その前にもう少しだけ、「原子核＋内殻電子の有効ポテンシャル」$V(\boldsymbol{r})$がどのようなものなのか考える必要がありそうです。

特に、tight-binding近似（和名：強束縛近似）はその名前の通り固体中の（内殻電子じゃない）電子を、「原子核に強く束縛された状態」つまり「ほぼ孤立原子中の電子」として扱う考え方です。

**というわけで、本章では水素様原子中から一歩進んで、「原子核＋多数の電子」の状態つまり「多電子原子」の状態の考え方について整理していきます。**

とはいえ、これも色々と歴史が長く研究レベルの内容まで踏み込むのが厳しかったので、具体的に解けるもっとも簡単な多電子原子、つまりHe原子（原子核＋2電子）の状態を中心に扱い、多電子原子の場合の考え方について雰囲気をつかめればOK、という方針で書いていくことにします。



# He原子

まずはもっとも簡単な多電子原子、He原子のシュレディンガー方程式を題材に、多電子原子の波動関数を理解していきます。例によって雑なイラストで表すと、

![](/images/tb/he.png)

こんな感じです。

[水素様原子の場合](https://zenn.dev/ponzumai/articles/tight-binding-model-hydrogen-atom)は、重心運動と相対運動に分離することにより、原子核位置を原点とした電子の相対座標のシュレディンガー方程式を導出しました。
しかし今回は3体の座標を考える必要があります。この場合にきれいに運動を分離できるかどうか私は知りません。

3体以上が絡む運動を考える場合、思い切って原子核が静止していることにして問題を解いていくことにします。これを **ボルンーオッペンハイマー近似（ Born–Oppenheimer approximation）** と言います。大層なことを書きましたが、とにかく原子核の運動は考えないことにして、原子核位置を原点とし、2電子の運動のみを考えていくことにします。

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

まずは前章の復習も兼ねて、2体相互作用項を無視するところから始めてみます。

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

## 平均場近似とSelf-consistent法

上記電子間相互作用を無視した解を求めましたが、実際には電子はそれぞれもう一方の電子からの斥力ポテンシャルを感じているはずです。
特に基底状態の場合に対して、その効果を取り入れた近似的な波動関数の形を求める方法が **平均場近似** （特に今回はHartree近似）と **self-consistent法**（自己無撞着法）です。これは後に出てくる、一般の多電子状態に適用するHartree-Fock近似の最も簡単な場合です。

<!--参考：http://cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/atom.html-->

近似的な波動関数を求めるために、量子力学の変分法を用います。
すなわち、試行関数を用意して、ハミルトニアン演算子の期待値、つまりエネルギー期待値を最小にするという条件を満たす関数形を求めます。

具体的には、基底状態の波動関数の解の形を何らかの軌道関数$\varphi(\boldsymbol{r})$にスピン関数$\alpha(\sigma), \beta(\sigma)$をそれぞれかけ合わせた関数からなるスレーター行列式と仮定します。これを試行関数と呼びます。

特に、今回はスピン上向き・下向きで同じ状態を2つまでとれるので、軌道部分は同じ関数形を仮定します。

すなわち試行関数を

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

と設定し、

ハミルトニアン演算子

$$
\mathcal{H} = \frac{-\hbar^2}{2m}\nabla_1{}^2 + \frac{-\hbar^2}{2m}\nabla_2{}^2 
+ \frac{-Ze^2}{4\pi\varepsilon_0r_1}
+ \frac{-Ze^2}{4\pi\varepsilon_0r_2}
+ \frac{e^2}{4\pi\varepsilon_0r_{12}}
$$

に対するエネルギー期待値を最小にするような関数$\varphi$の条件を導きます。

具体的にエネルギー期待値を計算していきましょう。まずはハミルトニアンにスピンへの作用がないので、スピン部分だけを独立で積分します。途中、$\varphi(\boldsymbol{r_1}), \alpha(\sigma_1)$を$\varphi(1)\alpha(1)$等と略記して、

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

続いて軌道部分はハミルトニアン演算子の中身を具体的に書いて展開すると、

$$
\begin{align*}
&\int \varphi^*(1)\varphi^*(2)\mathcal{H}\varphi(1)\varphi(2)d\boldsymbol{r}_1d\boldsymbol{r}_2\\
&=
    2\int\varphi^*(\boldsymbol{r})\left( \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-Ze^2}{4\pi\varepsilon_0r} \right)
    \varphi(\boldsymbol{r})d\boldsymbol{r}\\
&\>\>\>\>+
    \int\int \varphi^*(\boldsymbol{r}_1)\varphi^*(\boldsymbol{r}_2)\frac{e^2}{4\pi\varepsilon_0\left|\boldsymbol{r}_{1}-\boldsymbol{r}_2\right|}
    \varphi(\boldsymbol{r}_1)\varphi(\boldsymbol{r}_2)d\boldsymbol{r}_1d\boldsymbol{r}_2\\
&\equiv
    2\int\varphi^*(\boldsymbol{r})\hat{H}_1
    \varphi(\boldsymbol{r})d\boldsymbol{r}
    +
    \int\int \varphi^*(\boldsymbol{r}_1)\varphi^*(\boldsymbol{r}_2)
    V(\boldsymbol{r}_1,\boldsymbol{r}_2)
    \varphi(\boldsymbol{r}_1)\varphi(\boldsymbol{r}_2)d\boldsymbol{r}_1d\boldsymbol{r}_2
\end{align*}
$$

となります。ここでハミルトニアンの一体部分を$\hat{H}_1$、二体部分を$V(\boldsymbol{r}_1,\boldsymbol{r}_2)$と置きました。

ここに、軌道関数の規格化条件

$$
\int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}=1
$$

を束縛条件とし、未定乗数$\varepsilon$を用いて、次の汎関数を最小にする関数$\varphi$を見つければよいことになります。軌道関数$\varphi$と複素共役$\varphi^*$を異なる関数と考えて汎関数を$I[\varphi^*,\varphi]$と置くと、

$$
I\left[\varphi^*,\varphi\right] = 
    2\int\varphi^*(\boldsymbol{r})\hat{H}
    \varphi(\boldsymbol{r})d\boldsymbol{r}
    +
    \int \int\varphi^*(\boldsymbol{r}_1)\varphi^*(\boldsymbol{r}_2)
    V(\boldsymbol{r}_1,\boldsymbol{r}_2)
    \varphi(\boldsymbol{r}_1)\varphi(\boldsymbol{r}_2)d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
    -\varepsilon
    \left[
        \int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}-1  
    \right]
$$

となり、この変分は

$$
\begin{align*}
\delta I &= 
    I\left[\varphi^* + \delta\varphi^*, \varphi + \delta\varphi\right]
    - I\left[\varphi^*,\varphi\right]\\

&=
    2\int\left(
        \varphi^*(\boldsymbol{r})+\delta\varphi^*(\boldsymbol{r})
    \right)
    \hat{H}
    \left(
        \varphi(\boldsymbol{r})+\delta\varphi(\boldsymbol{r})
    \right)d\boldsymbol{r}\\
&\>\>\>\>+
    \int \left(
        \varphi^*(\boldsymbol{r}_1)+\delta\varphi^*
    \right)
    \left(
        \varphi^*(\boldsymbol{r}_2)+\delta\varphi^*
    \right)
    V(\boldsymbol{r}_1,\boldsymbol{r}_2)
    \left(
        \varphi(\boldsymbol{r}_1)+\delta\varphi
    \right)
    \left(
        \varphi(\boldsymbol{r}_2)+\delta\varphi
    \right)
    d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
    &\>\>\>\>-\varepsilon
    \left[
        \int
        \left(
            \varphi^*(\boldsymbol{r})+\delta\varphi^*
         \right)
        \left(
            \varphi(\boldsymbol{r})+\delta\varphi
        \right)
    d\boldsymbol{r}-1  
    \right]\\
&\>\>\>\>-
2\int\varphi^*(\boldsymbol{r})\hat{H}
    \varphi(\boldsymbol{r})d\boldsymbol{r}
    +
    \int \varphi^*(\boldsymbol{r}_1)\varphi^*(\boldsymbol{r}_2)
    V(\boldsymbol{r}_1,\boldsymbol{r}_2)
    \varphi(\boldsymbol{r}_1)\varphi(\boldsymbol{r}_2)d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
    &\>\>\>\>-\varepsilon
    \left[
        \int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}-1  
    \right]\\
&\>\>\>\>-I[\varphi^*,\varphi]

\end{align*}
$$
を考えれば良いことになります。
ここからかっこの中を展開して、微小変化の二次以上を含む項を落として整理すると、


$$
\begin{align*}

&\simeq
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
    d\boldsymbol{r}\\
&+
    \int 
      \delta\varphi^*(\boldsymbol{r}_1)
        \varphi^*(\boldsymbol{r}_2)
    V(\boldsymbol{r}_1,\boldsymbol{r}_2)
        \varphi(\boldsymbol{r}_1)
        \varphi(\boldsymbol{r}_2)
    d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
&+
 \int 
      \varphi^*(\boldsymbol{r}_1)
        \delta\varphi^*(\boldsymbol{r}_2)
    \mathcal{H}
        \varphi(\boldsymbol{r}_1)
        \varphi(\boldsymbol{r}_2)
    d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
&+
 \int 
    \varphi^*(\boldsymbol{r}_1)
        \varphi^*(\boldsymbol{r}_2)
    \mathcal{H}
        \delta\varphi(\boldsymbol{r}_1)
        \varphi(\boldsymbol{r}_2)
    d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
&+
 \int 
      \varphi^*(\boldsymbol{r}_1)
    \varphi^*(\boldsymbol{r}_2)
    \mathcal{H}
        \varphi(\boldsymbol{r}_1)
        \delta\varphi(\boldsymbol{r}_2)
    d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
    &\>\>\>\>-\varepsilon\int
        \delta\varphi^* \varphi(\boldsymbol{r}) d\boldsymbol{r}
    -
    \varepsilon\int
            \varphi^*(\boldsymbol{r})\delta\varphi d\boldsymbol{r} \\

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
    d\boldsymbol{r}\\
 &+
    2\int 
      \delta\varphi^*(\boldsymbol{r}_1)
        \varphi^*(\boldsymbol{r}_2)
    V(\boldsymbol{r}_1,\boldsymbol{r}_2)
        \varphi(\boldsymbol{r}_1)
        \varphi(\boldsymbol{r}_2)
    d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
&+
 2\int 
    \varphi^*(\boldsymbol{r}_1)
    \varphi^*(\boldsymbol{r}_2)
    V(\boldsymbol{r}_1,\boldsymbol{r}_2)
        \delta\varphi(\boldsymbol{r}_1)
        \varphi(\boldsymbol{r}_2)
    d\boldsymbol{r}_1d\boldsymbol{r}_2 \\
    &\>\>\>\>-\varepsilon\int
        \delta\varphi^* \varphi(\boldsymbol{r}) d\boldsymbol{r}
    -
    \varepsilon\int
            \varphi^*(\boldsymbol{r})\delta\varphi d\boldsymbol{r} \\
&=
    2\int\delta\varphi^*(\boldsymbol{r})
    \left[
        \hat{H}
            \varphi(\boldsymbol{r})
        +
        \left(
        \int
            \varphi^*(\boldsymbol{r}')
            V(\boldsymbol{r},\boldsymbol{r}')
            \varphi(\boldsymbol{r}')  d\boldsymbol{r}'
        \right)
        \varphi(\boldsymbol{r})

        -\frac{\varepsilon}{2}\varphi(\boldsymbol{r})
    \right] d\boldsymbol{r}\\
&\>\>\>\>+
    2\int
    \left[
        \varphi^*(\boldsymbol{r})
        \hat{H}
        +
        \varphi^*(\boldsymbol{r})
        \left(
        \int
            \varphi^*(\boldsymbol{r}')
            V(\boldsymbol{r},\boldsymbol{r}')
            \varphi(\boldsymbol{r}')  d\boldsymbol{r}'
        \right)
        -\frac{\varepsilon}{2}\varphi^*(\boldsymbol{r})
    \right]
    \delta\varphi(\boldsymbol{r}) d\boldsymbol{r}
\end{align*}
$$

となります。ここで最後の式から、任意の微小変化$\delta\varphi^*, \delta\varphi$に対して積分がゼロになる事から、$[\>\>\>]$の中の関数がゼロとなる必要があり、等式が二つ出てきますが、
エルミート演算子の性質を用いて

$$
\left\{\varphi^*\hat{H}\right\}^\dagger
=
\hat{H}^\dagger\varphi
=
\hat{H}\varphi
$$

となることから、これら二つの方程式を同値とできて、
$\varepsilon/2 = \varepsilon$と置きなおすと、最終的に近似関数が満たすべき一つの方程式

$$
\hat{H} \varphi(\boldsymbol{r})
+
\left(
        \int
            \varphi^*(\boldsymbol{r}')
            V(\boldsymbol{r},\boldsymbol{r}')
            \varphi(\boldsymbol{r}')  d\boldsymbol{r}'
        \right)
        \varphi(\boldsymbol{r})
        -\varepsilon\varphi(\boldsymbol{r})

        =0\\
\Rightarrow
    \hat{H}
            \varphi(\boldsymbol{r})
        +
        \left(
        \int
            \varphi^*(\boldsymbol{r}')
            V(\boldsymbol{r},\boldsymbol{r}')
            \varphi(\boldsymbol{r}')  d\boldsymbol{r}'
        \right)
        \varphi(\boldsymbol{r})
        =\varepsilon\varphi(\boldsymbol{r})
$$

を導くことができました。

### 平均場近似・Hartreeの近似

上式は一見、普通のシュレーディンガー方程式のような形をしていますが、よく見ると求めるべき$\varphi$が演算子部分に入ってしまっています。ここで積分部分を見てみると、

$$
\int
\varphi^*(\boldsymbol{r}')
V(\boldsymbol{r},\boldsymbol{r}')
\varphi(\boldsymbol{r}')  d\boldsymbol{r}'
=
\int
V(\boldsymbol{r},\boldsymbol{r}')
\left|\varphi(\boldsymbol{r}') \right|^2 d\boldsymbol{r}'
$$

これは着目している電子から見て、もう一方の電子から受けるクーロン相互作用に、その電子の確率密度がかかった形になっています。

そこでその部分を先に積分してしまい、

$$
\int
V(\boldsymbol{r},\boldsymbol{r}')
\left|\varphi(\boldsymbol{r}') \right|^2 d\boldsymbol{r}'
\equiv
    V(\boldsymbol{r})
$$

とすると、この部分は「着目している電子がもう一方の電子から感じるポテンシャル」の形になっています。

これは物理的には、元々は2電子間の相互作用を表していた部分を、「一つの電子は雲のように空間中に分布しており、そのような「電荷雲」からのポテンシャルを感じながらもう一つの電子が運動している」
という描像で近似していることに対応しています。

具体的には、空間に雲のように広がった電荷分布のようなものを考え、各点の電荷密度が
$$
-e|\varphi(\boldsymbol{r}')|^2
$$
であり、全空間を集めると合計の電荷が$-e$になるような電荷雲を考えると、そのような電荷分布の微小体積$d\boldsymbol{r}$部分から位置$\boldsymbol{r}$の電子が受けるポテンシャルは

$$
\frac{e^2}{4\pi\epsilon_0|\boldsymbol{r}-\boldsymbol{r}'|}|\varphi(\boldsymbol{r}')|^2d\boldsymbol{r'}
$$

となり、これを全空間について積分することで、電荷雲から電子が受けるポテンシャルを考えることができます。

### self-consistent 法

上記のような考え方のもと、

$$
\hat{H}
            \varphi(\boldsymbol{r})
        +
        V(\boldsymbol{r})
        \varphi(\boldsymbol{r})
        =\varepsilon\varphi(\boldsymbol{r})
$$

という方程式の解を見つけることができればそれを用いてHe原子（2電子原子）の近似的な波動関数を作ることができます。ただし依然として$V(\boldsymbol{r})$を計算するためには関数$\varphi$が必要です。

このような演算子の中に求めるべき関数が含まれている方程式を「self-consistent方程式」とか「自己無撞着方程式」とか呼び、数値的な解法をself-consistent法などと呼びます。具体的には、以下のような数値計算を行うプログラムを書くことになります。

1. 初期値として適当な関数$\varphi_0(\boldsymbol{r})$を用意する
2.  $\varphi_0(\boldsymbol{r})$を用いて積分$\int\varphi^* (\boldsymbol{r}_2)\mathcal{H}_2\varphi(\boldsymbol{r}_2)  d\boldsymbol{r}_2$を計算し、$V(\boldsymbol{r})$を求める
3.  2.で求めた関数$V(\boldsymbol{r})$を用いた微分方程式$\left(    \hat{H}_1 + V(\boldsymbol{r})\right)\varphi(\boldsymbol{r})=\varepsilon\varphi(\boldsymbol{r})$を解いて固有関数$\varphi_1$と固有値$\varepsilon_1$を求める
4.  得られた$\varphi_1$を用いて2.,3.を行い、固有関数$\varphi_2$と固有値$\varepsilon_2$を得る
5.  上記を、固有関数$\varphi_n$と一つ前の固有関数$\varphi_{n-1}$の差が十分小さくなる、または固有値$\varepsilon_n$と一つ前の固有値$\varepsilon_{n-1}$の差が十分小さくなる：$|\varepsilon_n-\varepsilon_{n-1}|<\epsilon$まで続けて、そこで得られた固有関数$\varphi_n$を近似的な波動関数とする。

### ポテンシャルの球対称化

この際、基底状態の一体の固有状態を球対称な関数と仮定すれば、有効ポテンシャル

$$
\int
\frac{e^2}{4\pi\epsilon_0|\boldsymbol{r}-\boldsymbol{r}'|}|
\left|\varphi(|\boldsymbol{r}'|) \right|^2 d\boldsymbol{r}'
\equiv
    V(\boldsymbol{r})
$$

も$\boldsymbol{r}$の角度によらない球対称な関数

$$
V(|\boldsymbol{r}|) = V(r)
$$

となるはずです（具体的な計算は後ほど示します）。


したがって一体のシュレーディンガー方程式のポテンシャル部分が球対称関数となり、水素原子の章で行ったように動径部分の微分方程式と角度部分の微分方程式に分解することができるようになります。具体的には、

$$
\hat{H}
            \varphi(\boldsymbol{r})
        +
        V(r)
        \varphi(\boldsymbol{r})
        =\varepsilon\varphi(\boldsymbol{r})
$$

として

$$
\hat{H}=\left( \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-Ze^2}{4\pi\varepsilon_0r} \right)
$$

をばらして書き直すと、

$$
\left\{
    \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-Ze^2}{4\pi\varepsilon_0r}
    + V(r)
\right\}
\varphi(\boldsymbol{r}) = \varepsilon\varphi(\boldsymbol{r})
$$

となります。これは水素様原子中の電子のシュレーディンガー方程式

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


$$
-\frac{\hbar^2}{2m} \left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        \right)
        +
        \left(
        -\frac{Ze}{4\pi\epsilon _0r}
        + 
        \int
        \frac{e^2}{4\pi\epsilon_0|\boldsymbol{r}-\boldsymbol{r}'|}|
        \left|\varphi(|\boldsymbol{r}'|) \right|^2 d\boldsymbol{r}'
        \right) R = \epsilon R
$$

を自己無撞着に解けば良いことになりました。

### いったんまとめ

「はじめに」で、「「原子核＋内殻電子の有効ポテンシャル」$V(\boldsymbol{r})$がどのようなものなのか」と書きました。本節ではHe原子の例を扱ったので内殻電子ではありませんが、一般の多電子原子の場合でも同様に考えると **原子核からの引力クーロン相互作用と、他の電子からの斥力クーロン相互作用を平均したもの** という形に近似できそうだという見通しが得られました。

またこの方程式の解はポテンシャル項
$V(r) = \int\frac{e^2}{4\pi\epsilon_0|\boldsymbol{r}-\boldsymbol{r}'|}|\varphi(|\boldsymbol{r}'|)|^2 d\boldsymbol{r}'$
が加わっているため水素原子の場合とは異なる関数形になりますが、一体のポテンシャル部分は「核からの引力クーロンポテンシャル＋雲のように広がった球対称の斥力クーロンポテンシャル」という形であり、「電子雲からの斥力で弱まった引力クーロンポテンシャル」のようなモノになりそうです。

したがって、解となる波動関数も大体が水素原子の場合に求めた固有状態、「1s状態」とか「3d状態」のような状態に近いものになりそうです。

最後に、少し違う方法で上記の予想を確かめてみてHe原子の項を終わりにします。

### おまけ：遮蔽されたクーロンポテンシャル

先ほどと同様に基底状態の反対称化された波動関数を

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

と仮定し、かつ軌道関数の試行関数として,
原子核の電荷が$+Ze$の水素様原子中の電子の$1s$状態

$$
\varphi_{1s} = \varphi _{100} = R_{10}Y_0^0,\\
R_{10}(r) = \left( \frac{Z}{a_0}  \right)^{3/2}2e^{-Zr/a_0} \equiv Ae^{-Zr/a_0}\\
Y_0^0 = \frac{1}{\sqrt{4\pi}} 
$$

を使います。（動径関数の係数を簡略化のため$A$と置きました）


「原子核からのクーロンポテンシャルが弱まった状態」を予想しているので、原子核の電荷に対応する$Z$を変分法で最適化するパラメータと考え、変分計算を具体的に行っていきます。

ここで注意点としては、ハミルトニアンはHe原子のハミルトニアン、つまり$Z=2$として考える必要があります。

$$
\mathcal{H} = \frac{-\hbar^2}{2m}\nabla_1{}^2 + \frac{-\hbar^2}{2m}\nabla_2{}^2 
+ \frac{-2e^2}{4\pi\varepsilon_0r_1}
+ \frac{-2e^2}{4\pi\varepsilon_0r_2}
+ \frac{e^2}{4\pi\varepsilon_0r_{12}}
$$


1s状態は原点からの距離$|\boldsymbol{r}|=r$のみの関数なので以降$\varphi_{1s}(r)$と書き、エネルギー期待値$\langle E \rangle$は、スピン部分を先ほどと同じく先に積分して

$$
\begin{align*}
\langle E \rangle &= \int \varphi_{1s}^*(1)\varphi_{1s}^*(2)\mathcal{H}\varphi_{1s}(1)\varphi_{1s}(2)d\boldsymbol{r}_1d\boldsymbol{r}_2\\
&=
    2\int\varphi_{1s}^*(r)\left( \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-2e^2}{4\pi\varepsilon_0r} \right)
    \varphi_{1s}(r)d\boldsymbol{r}\\
&\>\>\>\>+
    \int\int \varphi_{1s}^*(r)\varphi_{1s}^*(r')\frac{e^2}{4\pi\varepsilon_0\left|\boldsymbol{r}-\boldsymbol{r}'\right|}
    \varphi_{1s}(r)\varphi_{1s}(r')d\boldsymbol{r}d\boldsymbol{r}'
\end{align*}
$$

初めに簡単な一体部分から考えていきます。^[水素様原子の一体のハミルトニアンに対して水素様原子の固有関数を考えているので、一見単純に固有値を使えばいいように思えますが（というか私が最初そのように勘違いして計算して答えが合わなかったのですが）、ハミルトニアンは$Z=2$、試行関数は$Z$を変数として扱っているので実は固有関数にならず、ちゃんと計算しないといけないのでした。]

はじめにポテンシャル部分について、$\varphi_{1s}$の関数を入れて角度部分を先に積分します。

$$
\begin{align*}

2\int\varphi_{1s}^*(r)\left(  \frac{-2e^2}{4\pi\varepsilon_0r} \right)
    \varphi_{1s}(r)d\boldsymbol{r}

&=

2\times\frac{-2e^22^2}{4\pi(4\pi\varepsilon_0)}
\left(\frac{Z}{a_0}\right)^3\int
    \frac{1}{r} 
    e^{-2Zr/a_0}
    d\boldsymbol{r}\\

&=
\frac{-2^4e^2}{4\pi(4\pi\varepsilon_0)}
\left(\frac{Z}{a_0}\right)^3
\int
    \frac{1}{r} 
    e^{-2Zr/a_0}
    r^2dr
    \iint\sin\theta d\theta d \phi\\

&=
\frac{-2^4e^2}{4\pi\varepsilon_0}
\left(\frac{Z}{a_0}\right)^3
\int
    re^{-2Zr/a_0}
    dr
\end{align*}
$$

ここで$\int r e^{-2\alpha r}dr$について以下の関係を考えて、

$$
-\frac{1}{2\alpha}
    \left\{
        \left(r+\frac{1}{2\alpha}\right)
    e^{-2\alpha r} 
    \right\}' = re^{-2\alpha r}
$$

定積分を計算すると、

$$
\begin{align*}
\frac{-2^4e^2}{4\pi\varepsilon_0}
\left(\frac{Z}{a_0}\right)^3
\int
    re^{-2Zr/a_0}
    dr
&=
\frac{-2^4e^2}{4\pi\varepsilon_0}
\left(\frac{Z}{a_0}\right)^3
    \frac{-1}{2}\left(\frac{Z}{a_0}\right)^{-1}
    \left[
    \left\{
        r + \frac{1}{2}\left(\frac{Z}{a_0}\right)^{-1}
    \right\}
    e^{-2Zr/a_0}   
    \right]_0^\infty\\

&=
\frac{2^3e^2}{4\pi\varepsilon_0}
\left(\frac{Z}{a_0}\right)^3
    \left(\frac{Z}{a_0}\right)^{-1}
    \left\{
        -\frac{1}{2}\left(\frac{Z}{a_0}\right)^{-1}
        \right\}\\

&=
    \frac{-4e^2}{4\pi\varepsilon_0}
    \frac{Z}{a_0}

\end{align*}
$$

続いて運動エネルギー部分は、極座標表示のラプラシアンの$e^{-Zr/a_0}$への作用

$$
\begin{align*}
\left(
    \frac{d^2}{dr^2} + \frac{2}{r}\frac{d}{dr} + 
    \frac{1}{r}\Lambda(\theta,\phi)
\right)
e^{-Zr/a_0}

&=
    
    \left(
        \frac{d^2}{dr^2} + \frac{2}{r}\frac{d}{dr} 
    \right)
    e^{-Zr/a_0}\\

&=
    \left\{
    \left(\frac{Z}{a_0}\right)^2
    +
    \frac{2}{r}\frac{Z}{a_0}
    \right\}e^{-Zr/a_0}
\end{align*}
$$



を用いて

$$
\begin{align*}

&2\int\varphi_{1s}^*(r)\left( 
    \frac{-\hbar^2}{2m}\nabla{}^2 
    \right)
    \varphi_{1s}(r)d\boldsymbol{r}\\


&=
    2 \frac{-\hbar^2}{2m}
    2^2\left(\frac{Z}{a_0}\right)^3
    \frac{1}{4\pi}
\int 
    \left\{
    \left(\frac{Z}{a_0}\right)^2
    -
    \frac{2}{r}\frac{Z}{a_0}
    \right\}e^{-2Zr/a_0}
    r^2dr
    \iint\sin\theta d\theta d\phi\\

&=
    \frac{-2^3\hbar^2}{2m}
    \left(\frac{Z}{a_0}\right)^3
    \int
    \left\{
    \left(\frac{Z}{a_0}\right)^2r^2
    -
    2r\frac{Z}{a_0}
    \right\}e^{-2Zr/a_0}
    dr
\end{align*}
$$

$$
-\frac{1}{2\alpha}
    \left\{
        \left(r+\frac{1}{2\alpha}\right)
    e^{-2\alpha r} 
    \right\}' = re^{-2\alpha r},\\
-\frac{1}{2\alpha}\left\{
    \left(
    r^2
    +
    \frac{1}{\alpha}r
    +
     \frac{1}{2\alpha^2} 
     \right)e^{-2\alpha r}
    \right\}' = r^2e^{-2\alpha r}

$$


$$
\begin{align*}
 &\frac{-2^3\hbar^2}{2m}
    \left(\frac{Z}{a_0}\right)^3
    \int
    \left\{
    \left(\frac{Z}{a_0}\right)^2r^2
    -
    2r\frac{Z}{a_0}
    \right\}e^{-2Zr/a_0}
    dr\\
&=
 \frac{-2^3\hbar^2}{2m}
    \left(\frac{Z}{a_0}\right)^3
    \left(
        \frac{-1}{2}\left(\frac{Z}{a_0}\right)\left[
            \left\{
                r^2 +\left( \frac{Z}{a_0}\right)^{-1}r + 
                \frac{1}{2}\left( \frac{Z}{a_0}\right)^{-2}
            \right\}e^{-2Zr/a_0}
            \right]_0^{\infty}\right.\\            
&\>\>\>\>+
\left.
        
        \left[
            \left\{
                r + \frac{1}{2}\left( \frac{Z}{a_0}\right)^{-1}
                \right\}e^{-2Zr/a_0}
        \right]_0^{\infty}
        

    \right)\\

&=
    \frac{-2^3\hbar^2}{2m}
    \left(\frac{Z}{a_0}\right)^3
    \left[
        \frac{1}{4}\left(\frac{Z}{a_0}\right)^{-1}
        -
        \frac{1}{2}\left(\frac{Z}{a_0}\right)^{-1}
        \right]\\

&=
    \frac{-\hbar^2}{m}\left(\frac{Z}{a_0}\right)^2
\end{align*}
$$

$$
a_0 = \frac{4\pi\varepsilon_0\hbar^2}{me^2}
$$

$$
\frac{-\hbar^2}{m}\left(\frac{Z}{a_0}\right)^2
=
\frac{-e^2}{4\pi\varepsilon_0}\frac{Z^2}{a_0}
$$


となります。

$\langle E \rangle$の第2項

$$
\iint \varphi_{1s}^*(r)\varphi_{1s}^*(r')\frac{e^2}{4\pi\varepsilon_0\left|\boldsymbol{r}-\boldsymbol{r}'\right|}
    \varphi_{1s}(r)\varphi_{1s}(r')d\boldsymbol{r}d\boldsymbol{r}'
$$
について、先に$\boldsymbol{r}'$についての積分を行います。これは有効ポテンシャル$V(r)$に対応します^[この部分の計算は[こちらのノート](http://www.th.phys.titech.ac.jp/~muto/lectures/QMII10/QMII10_chap17.pdf)を参考にしました。感謝です。]。

積分では$\boldsymbol{r}$を固定して$\boldsymbol{r}'$についての積分を求めるので、以下のように
$\boldsymbol{r}$の方向を$z$軸、にとって極座標表示を行うと、
![](/images/tb/rpmr.png)

$$
|\boldsymbol{r}' - \boldsymbol{r}| = \sqrt{(r-r'\cos\theta)^2+(r'\sin\theta)^2}\\
=\sqrt{r^2-2rr'\cos\theta+r'^2}
$$

となるので、

$$
\begin{align*}
\int
\frac{1}{|\boldsymbol{r}-\boldsymbol{r}'|}|
\left|\varphi_{1s}(r') \right|^2 d\boldsymbol{r}'
&=
\iiint\frac{1}{\sqrt{r^2+2rr'\cos\theta'+r'^2}}
\frac{A^2}{4\pi}e^{-2Zr'/a_0}r'^2\sin\theta' dr'd\theta' d\phi'
\end{align*}
$$

となる。非積分関数は$\phi$によらないので独立で積分でき、

$$
\int_0^{2\pi} d\phi' = 2\pi
$$

です。

また$\theta'$に関する積分は、$t = \cos\theta'$と置換して$dt = -\sin\theta' d\theta'$より

$$
\begin{align*}
\int_0^{\pi}\frac{1}{\sqrt{r^2-2rr'\cos\theta'+r'^2}}
\sin\theta' d\theta' 
&= -\int_{-1}^{1}\frac{1}{\sqrt{r^2+r'^2-2rr't}}dt\\
&=
    \frac{1}{rr'}\left[
        \sqrt{r^2+r'^2-2rr't}
        \right]_{-1}^{1}\\
&=
    \frac{1}{rr'}\left\{
        \sqrt{(r-r')^2} - \sqrt{(r+r')^2}
        \right\}\\
&=
    \frac{1}{rr'}\left\{
        |r-r'| - r+r'
        \right\}\\
&=
\left\{
\begin{array}{ll}
\frac{2}{r} & r'<r \\
\\
\frac{2}{r'} & r'>r
\end{array}
\right.

\end{align*}
$$

となります。$\phi'$の積分で出てきた係数をかけて最後に残った$r'$に関する積分を書くと、積分範囲を$r>r'とr<r'$の部分に分けて

$$
\begin{align*}
\int
\frac{1}{|\boldsymbol{r}-\boldsymbol{r}'|}|
\left|\varphi_{1s}(r') \right|^2 d\boldsymbol{r}'

&=
\int_0^r
\frac{A^2}{r}e^{-2Zr'/a_0}r'^2dr'

+
\int_r^{\infty}
\frac{A^2}{r'}e^{-2Zr'/a_0}r'^2dr'\\

&=
\frac{A^2}{r}\int_0^r
e^{-2Zr'/a_0}r'^2dr'

+
A^2\int_r^{\infty}
e^{-2Zr'/a_0}r'dr'
\end{align*}
$$

この積分は

$$
-\frac{1}{\alpha}
    \left\{
        \left(r+\frac{1}{\alpha}\right)
    e^{-\alpha r} 
    \right\}' = re^{-\alpha r},\\
-\frac{1}{\alpha}\left\{
    \left(
    r^2
    +
    \frac{2}{\alpha}r
    +
     \frac{2}{\alpha^2} 
     \right)e^{-\alpha r}
    \right\}' = r^2e^{-\alpha r}
$$

を考えて定積分の計算をして
$$
A = 2\left(\frac{Z}{a_0}\right)^{3/2}
$$

を戻すと、最終的に


$$
\begin{align*}
\int
\frac{e^2}{4\pi\varepsilon_0|\boldsymbol{r}-\boldsymbol{r}'|}|
\left|\varphi_{1s}(r') \right|^2 d\boldsymbol{r}'

=
\frac{e^2}{4\pi\varepsilon_0r}
    \left\{
        1-\left(1+\frac{Z}{a_0}r\right)e^{-2Zr/a_0}
        \right\}
\end{align*}
$$

となります。今度はこの関数を用いて$r$に関する積分

$$
\int
\frac{e^2}{4\pi\varepsilon_0r}
    \left\{
        1-\left(1+\frac{Z}{a_0}r\right)e^{-2Zr/a_0}
        \right\}
    \left|\varphi_{1s}(r) \right|^2d\boldsymbol{r}
$$

を計算すれば$\langle E \rangle$を求めることができますが、ここで少し立ち止まって有効ポテンシャルの積分の結果を眺めてみると、確かに「大体クーロンポテンシャル」みたいな形をしており、特に$r\rightarrow\infty$では$\{\>\>\>\>\}$の中の第2項がゼロになり$Z=1$のクーロンポテンシャルになっていることがわかります。He原子の原子核の電荷は$Z=2$なので、$+2e$の正電荷を原点近傍に分布している$1s$軌道の電荷$-e$の「電荷雲」が遮蔽しているというイメージがより明確になりました。

さて、話を戻して$\boldsymbol{r}$に関する積分を行います。

非積分関数が球対称なので$\theta, \phi$に関する積分

$$
\iint \sin\theta d\theta d\phi = 4\pi
$$

を先に行い、最後に動径$r$の積分を考えると、先ほどと同じ形の積分が出てくるので同様に計算して、

$$
\begin{align*}
&\>\>\>\>\>\>4\left(
    \frac{Z}{a_0}
    \right)^3\int
\frac{e^2}{4\pi\varepsilon_0r}
    \left\{
        1-\left(1+\frac{Z}{a_0}r\right)e^{-2Zr/a_0}
        \right\}
    e^{-2Zr/a_0}r^2
    dr\\
&=
4\left(
    \frac{Z}{a_0}
    \right)^3
    \frac{e^2}{4\pi\varepsilon_0}
\int_0^\infty
\left(
    r e^{-2Zr/a_0}
    -
    r e^{-4Zr/a_0}
    -
    \frac{Z}{a_0}r^2 e^{-4Zr/a_0}
    \right) dr\\
&=
4\left(
    \frac{Z}{a_0}
    \right)^3
    \frac{e^2}{4\pi\varepsilon_0}
    \frac{5}{32}\left(
        \frac{Z}{a_0}
        \right)^{-2}\\
&=
\frac{e^2}{4\pi\varepsilon_0}\frac{5}{8}\frac{Z}{a_0}
\end{align*}
$$


と、気持ちよいくらいきれいな結果が得られました。

最後に随分前に計算した第1項

$$
2\varepsilon_1
= -\frac{Z^2e^2}{4\pi\varepsilon_0 a_0}
$$

と合わせて、

$$
\begin{align*}
\langle E \rangle
&=
-\frac{Z^2e^2}{4\pi\varepsilon_0 a_0}
+
\frac{e^2}{4\pi\varepsilon_0}\frac{5}{8}\frac{Z}{a_0}\\

&=
-\frac{e^2}{4\pi\varepsilon_0 a_0}\left(
    Z^2 - \frac{5}{8}Z
    \right)
\end{align*}
$$

となりました。この期待値を最小にする$Z$は、

$$
\frac{d\langle E \rangle}{dZ} = -\frac{e^2}{4\pi\varepsilon_0 a_0}\left(
    2Z - \frac{5}{8}
    \right)
= 0
$$

より、
$$
Z = \frac{5}{16}
$$

とすれば良いことになります。






# 一般の多電子原子
## Hartree近似

## Hartree-Fock近似

## 多電子原子内の電子の波動関数と平均場ポテンシャル


この際、微分方程式

$$
\hat{H}
            \varphi(\boldsymbol{r})
        +
        V(\boldsymbol{r})
        \varphi(\boldsymbol{r})
        =\varepsilon\varphi(\boldsymbol{r})
$$

は$V(\boldsymbol{r})$が球対称な形になっておりません。なので水素原子の章で行ったように、動径部分と角度部分に分解することができません。

しかし$V(\boldsymbol{r})$を角度部分で平均化することで、近似的に
球対称ポテンシャルと考え、水素原子の章で行ったように動径部分の微分方程式と角度部分の微分方程式に分解することができるようになります。

具体的には、位置$\boldsymbol{r}$のポテンシャル$V(\boldsymbol{r})$について、角度部分について積分した後に表面積$4\pi r^2$で割ることで、

$$
V(r) = \frac{1}{4\pi r^2}\int V(\boldsymbol{r})r^2\sin\theta d\theta d\phi\\
= \frac{1}{4\pi}\int V(\boldsymbol{r})\sin\theta d\theta d\phi
$$

と近似することで、めでたく球対称な微分方程式

$$
\hat{H}
            \varphi(\boldsymbol{r})
        +
        V(r)
        \varphi(\boldsymbol{r})
        =\varepsilon\varphi(\boldsymbol{r})
$$

にたどり着きました。