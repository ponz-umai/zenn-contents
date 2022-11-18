---
title: "多電子原子内の電子の波動関数と平均場ポテンシャル"
emoji: "🌌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
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
特に、tight-binding近似はその名前の通り固体中の電子を、「原子核に強く束縛された状態」つまり「ほぼ孤立原子中の電子」として扱う考え方です。

そこで固体物理に行く前の最後のステップとして、「原子核＋多数の電子」の状態つまり「多電子原子」の状態の考え方を理解する必要があります。

本章では具体的に解けるもっとも簡単な多電子原子、つまりHe原子（原子核＋2電子）の状態を中心に扱いながら、
多電子原子の場合の考え方について雰囲気を掴めることを目的に書いていこうと思います。


# He原子中の電子

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
\epsilon _n = -\frac{Z^2me^4}{(4\pi\epsilon_0 )^22\hbar^2} \frac{1}{n^2} , \>\>\>\> n=1,2,\cdots
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

その効果を一体の相互作用の形で近似的に取り入れ、上手いこと一体のシュレディンガー方程式に落とし込んで近似的な波動関数の形を求める方法が **平均場近似** （特に今回はHartree近似）と **self-consistent法**（自己無撞着法）です。これは後に出てくる、一般の多電子状態に適用するHartree-Fock近似の最も簡単な場合です。

<!--参考：http://cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/atom.html-->

### 変分法の超概要

近似的な波動関数を求めるために、量子力学の変分法を用います。
詳細はいつか記事にすることにして、適当な量子力学の教科書やWEB上のテキスト^[例えば、[このノート(pdf)](http://rokamoto.sakura.ne.jp/education/quantum/variation20200811B.pdf)や[このノート(pdf)](http://kurasawa.c.ooco.jp/qm_a.pdf)が比較的詳しいです。ただ、もっと初歩的な部分が省略されていたりするので、自分の理解のためにいつかちゃんとまとめたいと思っています。]を参照していただくことにして、概要としては

- 正しい波動関数はハミルトニアン演算子の期待値（つまりエネルギー期待値）を最小にするはずだという前提のもと
- 正しい関数に近そうな関数（試行関数）を仮定して、ハミルトニアンの期待値を最小にする条件を導く
- 上記条件をもとに近似的な波動関数の形を求める

というようなことを行います。また、波動関数は規格化されている必要があり、その条件をラグランジュの未定乗数法という形で取り入れます。

### He原子中への適用

具体的にHe原子中の電子に変分法を適用していきます。

相互作用を無視した場合でみたように、基底状態であればスピン上向き・下向きで同じ軌道関数をとれるので、He原子中の2つの電子は何らかの一体固有エネルギーが最低の、同じ状態を取っているものと期待されます。
そこで基底状態の波動関数の解の形を1つの軌道関数$\varphi(\boldsymbol{r})$にスピン関数$\alpha(\sigma), \beta(\sigma)$をそれぞれかけ合わせた関数からなるスレーター行列式と仮定し、これを試行関数とします。

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

と設定し、ハミルトニアン演算子

$$
\mathcal{H} = \frac{-\hbar^2}{2m}\nabla_1{}^2 + \frac{-\hbar^2}{2m}\nabla_2{}^2 
+ \frac{-Ze^2}{4\pi\varepsilon_0r_1}
+ \frac{-Ze^2}{4\pi\varepsilon_0r_2}
+ \frac{e^2}{4\pi\varepsilon_0r_{12}}
$$

に対するエネルギー期待値を最小にするような関数$\varphi$の条件を導きます。

具体的にエネルギー期待値を計算していきましょう。まずはハミルトニアンにスピンへの作用がないので、スピン部分だけを独立で積分します。
途中、$\varphi(\boldsymbol{r_1}), \alpha(\sigma_1)$を$\varphi(1)\alpha(1)$等と略記して、

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
    \int\varphi^*(\boldsymbol{r}_1)\left( \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-Ze^2}{4\pi\varepsilon_0r_1} \right)
    \varphi(\boldsymbol{r}_1)d\boldsymbol{r}_1
    
    \int\varphi^*(\boldsymbol{r}_2)
    \varphi(\boldsymbol{r}_2)d\boldsymbol{r}_2
    \\
&\>\>\>\>+

\int\varphi^*(\boldsymbol{r}_1)
    \varphi(\boldsymbol{r}_1)d\boldsymbol{r}_1
    
    \int\varphi^*(\boldsymbol{r}_2)
    \left( \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-Ze^2}{4\pi\varepsilon_0r_2} \right)
    \varphi(\boldsymbol{r}_2)d\boldsymbol{r}_2
    \\
&\>\>\>\>+


    \int\int \varphi^*(\boldsymbol{r}_1)\varphi^*(\boldsymbol{r}_2)\frac{e^2}{4\pi\varepsilon_0\left|\boldsymbol{r}_{1}-\boldsymbol{r}_2\right|}
    \varphi(\boldsymbol{r}_1)\varphi(\boldsymbol{r}_2)d\boldsymbol{r}_1d\boldsymbol{r}_2\\
&\equiv
    2\int\varphi^*(\boldsymbol{r})\hat{H}
    \varphi(\boldsymbol{r})d\boldsymbol{r}

    \int\varphi^*(\boldsymbol{r}')
    \varphi(\boldsymbol{r}')d\boldsymbol{r}'
    +
    \int\int \varphi^*(\boldsymbol{r})\varphi^*(\boldsymbol{r}')
    V(\boldsymbol{r},\boldsymbol{r}')
    \varphi(\boldsymbol{r})\varphi(\boldsymbol{r}')d\boldsymbol{r}d\boldsymbol{r}'
\end{align*}
$$

となります。ここで軌道部分の関数を同じ関数と仮定したので、$\boldsymbol{r}_1$と$\boldsymbol{r}_2$に関する積分が同じ形になります。そこで$\boldsymbol{r}_1$と$\boldsymbol{r}_2$という表記から$\boldsymbol{r}$と$\boldsymbol{r}'$という表記へ変更します。
また、ハミルトニアンの一体部分を$\hat{H}$、二体部分を$V(\boldsymbol{r},\boldsymbol{r}')$と置きました。

ここに、軌道関数の規格化条件

$$
\iint\left|\Phi(\boldsymbol{r}_1,\boldsymbol{r}_2)\right|^2d\boldsymbol{r}_1d\boldsymbol{r}_2=

\int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}\int\left|\varphi(\boldsymbol{r}')\right|^2d\boldsymbol{r}'=
1
$$

を束縛条件とし、未定乗数$E$を用いて、次の汎関数を最小にする関数$\varphi$を見つければよいことになります。

$$
I\left[\varphi\right] = 
    2\int\varphi^*(\boldsymbol{r})\hat{H}
    \varphi(\boldsymbol{r})d\boldsymbol{r}

    \int\varphi^*(\boldsymbol{r}')
    \varphi(\boldsymbol{r}')d\boldsymbol{r}'
    +
    \int\int \varphi^*(\boldsymbol{r})\varphi^*(\boldsymbol{r}')
    V(\boldsymbol{r},\boldsymbol{r}')
    \varphi(\boldsymbol{r})\varphi(\boldsymbol{r}')d\boldsymbol{r}d\boldsymbol{r}' \\
    -E
    \left[
        \int\left|\varphi(\boldsymbol{r})\right|^2d\boldsymbol{r}
        \int\left|\varphi(\boldsymbol{r}')\right|^2d\boldsymbol{r}'
        -1  
    \right]
$$

$\varphi$は複素関数なので$\varphi$と複素共役$\varphi^*$
のそれぞれに対する変分を考えればよく、$\varphi^*$についての変分を考えると、

$$
\begin{align*}
\delta I &= 
    I\left[\varphi^* + \delta\varphi^*, \varphi\right]
    - I\left[\varphi^*,\varphi\right]\\

&=
    2\int\left(
        \varphi^*(\boldsymbol{r})+\delta\varphi^*(\boldsymbol{r})
    \right)
    \hat{H}
    
        \varphi(\boldsymbol{r})
    d\boldsymbol{r}
    \int \varphi(\boldsymbol{r}')^*
    \varphi(\boldsymbol{r}')
    d\boldsymbol{r}'
    \\
&\>\>\>\>+
    \int \left(
        \varphi^*(\boldsymbol{r})+\delta\varphi^*
        (\boldsymbol{r})

    \right)
    \left(
        \varphi^*(\boldsymbol{r}')+\delta\varphi^*(\boldsymbol{r}')
    \right)
    V(\boldsymbol{r},\boldsymbol{r}')
    
        \varphi(\boldsymbol{r})
        \varphi(\boldsymbol{r}')
    
    d\boldsymbol{r}d\boldsymbol{r}' \\
    &\>\>\>\>-E
    \left[
        \int
        \left(
            \varphi^*(\boldsymbol{r})+\delta\varphi^*(\boldsymbol{r})
         \right)
         \varphi(\boldsymbol{r})
         d\boldsymbol{r}
         \int
         \left(
            \varphi^*(\boldsymbol{r}')+\delta\varphi^*(\boldsymbol{r}')
        \right)
        \varphi(\boldsymbol{r}')
    d\boldsymbol{r}'-1  
    \right]\\
&\>\>\>\>- I\left[\varphi^*,\varphi\right]
\end{align*}
$$

を考えれば良いことになります。
ここからかっこの中を展開して、微小変化の二次以上を含む項を落とし、また仮定した規格化条件より

$$
\int \left|\varphi(\boldsymbol{r}')\right|^2d\boldsymbol{r}' = 1
$$

とすると、変分は以下の形になります。


$$
\begin{align*}
\delta I &= 
    I\left[\varphi^* + \delta\varphi^*, \varphi\right]
    - I\left[\varphi^*,\varphi\right]\\

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

        -\left(
            E - \int\varphi^*(\boldsymbol{r}')
        \hat{H}
            \varphi(\boldsymbol{r}')
            d\boldsymbol{r}'
        \right)
            \varphi(\boldsymbol{r})
    \right] d\boldsymbol{r}\\

\end{align*}
$$

ここで任意の微小変化$\delta\varphi^*$に対して積分がゼロになる事から、$[\>\>\>]$の中の関数がゼロとなる必要があり、試行関数が満たすべき方程式を得ることができます。

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
        -\left(
            E - \int\varphi^*(\boldsymbol{r}')
        \hat{H}
            \varphi(\boldsymbol{r}')
            d\boldsymbol{r}'
        \right)
            \varphi(\boldsymbol{r})
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

ここで積分$\int \varphi^*(\boldsymbol{r}') \hat{H} \varphi(\boldsymbol{r}')  d\boldsymbol{r}'$
は定数となるので

$$
E - \int \varphi^*(\boldsymbol{r}') \hat{H} \varphi(\boldsymbol{r}')  d\boldsymbol{r}'
\equiv
\varepsilon
$$

と置くと、最終的に一体の波動関数が従う方程式は以下の固有値方程式の形となります。

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

なお、$\varphi$に関する変分からはこちらのエルミート共役となる方程式が得られますが、ハミルトニアンのエルミート性より同値な方程式となります。^[この「シュレディンガー方程式」の固有値$\varepsilon$は、$E-\int\varphi^*\hat{H}\varphi d\boldsymbol{r}$でした。この$E$は、二体の波動関数のシュレディンガー方程式$\mathcal{H}\Phi = E\Phi$の固有値に対応しています。（$\Phi$を一つの関数としてその変分を考え、同じように計算をすることで導けます）すなわち、一体の微分方程式の固有値$\varepsilon$は、「2電子の固有エネルギーから1電子分の固有エネルギーを差し引いたもの」の形をしていることがわかります。なお、全エネルギーを求めるためには$2\varepsilon \neq E$で、二重に数えている二体相互作用のポテンシャルを差し引く必要があります。つまり
$E = 2\varepsilon - \iint \varphi^*(\boldsymbol{r})\varphi^*(\boldsymbol{r}')V(\boldsymbol{r},\boldsymbol{r}') \varphi(\boldsymbol{r})\varphi(\boldsymbol{r}')d\boldsymbol{r} d\boldsymbol{r}'$。]


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
        =\varepsilon\varphi(\boldsymbol{r}),\\
\hat{H} = 
\frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-Ze^2}{4\pi\varepsilon_0r} ,\\

V(\boldsymbol{r}) = 
\int
V(\boldsymbol{r},\boldsymbol{r}')
\left|\varphi(\boldsymbol{r}') \right|^2 d\boldsymbol{r}'

$$

という方程式の解を見つけることができればそれを用いてHe原子（2電子原子）の近似的な波動関数を作ることができます。ただし依然として$V(\boldsymbol{r})$を計算するためには関数$\varphi$が必要です。

このような演算子の中に求めるべき関数が含まれている方程式を「self-consistent方程式」とか「自己無撞着方程式」とか呼び、数値的な解法をself-consistent法などと呼びます。具体的には、以下のような数値計算を行うプログラムを書くことになります。

1. 初期値として適当な関数$\varphi_0(\boldsymbol{r})$を用意する
2.  $\varphi_0(\boldsymbol{r})$を用いて積分$\int V(\boldsymbol{r},\boldsymbol{r}')\left|\varphi_0(\boldsymbol{r}') \right|^2 d\boldsymbol{r'}$を計算し、$V(\boldsymbol{r})$を求める
3.  2.で求めた関数$V(\boldsymbol{r})$を用いた微分方程式$\left(    \hat{H} + V(\boldsymbol{r})\right)\varphi(\boldsymbol{r})=\varepsilon\varphi(\boldsymbol{r})$を解いて固有関数$\varphi_1$と固有値$\varepsilon_1$を求める
4.  得られた$\varphi_1$を用いて2.,3.を行い、固有関数$\varphi_2$と固有値$\varepsilon_2$を得る
5.  上記を、固有関数$\varphi_n$と一つ前の固有関数$\varphi_{n-1}$の差が十分小さくなるまで続けて、そこで得られた固有関数$\varphi_n$を近似的な波動関数とする。^[固有値$\varepsilon_n$と一つ前の固有値$\varepsilon_{n-1}$の差が十分小さくなる：$|\varepsilon_n-\varepsilon_{n-1}|<\epsilon$とすることもある]

### 球対称な波動関数を仮定した場合

この際、基底状態の一体の固有状態を球対称な関数と仮定すれば、有効ポテンシャル

$$
\int
\frac{e^2}{4\pi\epsilon_0|\boldsymbol{r}-\boldsymbol{r}'|}|
\left|\varphi(|\boldsymbol{r}'|) \right|^2 d\boldsymbol{r}'
$$

も$\boldsymbol{r}$の角度によらない球対称な関数

$$
V(|\boldsymbol{r}|) = V(r)
$$

となるはずです（具体的な計算は後ほど示します）。


したがって一体のシュレーディンガー方程式のポテンシャル部分が球対称関数となり、水素原子の章で行ったように動径部分の微分方程式と角度部分の微分方程式に分解することができるようになります。

具体的には、求めるべき微分方程式は


$$
\left\{
    \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-Ze^2}{4\pi\varepsilon_0r}
    + V(r)
\right\}
\varphi(r) = \varepsilon\varphi(r)
$$

ですが、ラプラシアンの極座標表示

$$
\nabla^2 = 
\left(
    \frac{d^2}{dr^2} + \frac{2}{r}\frac{d}{dr} + 
    \frac{1}{r}\Lambda(\theta,\phi)
\right),\\

\hat{\Lambda} (\theta ,\phi) = \left\{
     \frac{1}{\sin \theta}  \frac{\partial }{\partial \theta } \left( \sin \theta \frac{\partial }{\partial \theta }  \right) 
     +
     \frac{1}{\sin ^2\theta }\frac{\partial ^2}{\partial \phi ^2}  
\right\} 
$$

を代入して、第3項は球対称な関数に対してはゼロになることから、

最終的に動径部分のみの1変数微分方程式

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

を自己無撞着に解けば良いことになります。

### いったんまとめ

「はじめに」で、「「原子核＋内殻電子の有効ポテンシャル」$V(\boldsymbol{r})$がどのようなものなのか」と書きました。本節ではHe原子の例を扱ったので内殻電子ではありませんが、一般の多電子原子の場合でも同様に考えると **①原子核からの引力クーロン相互作用ポテンシャル、②他の電子からの斥力クーロン相互作用を平均したポテンシャル** という形に近似できそうだという見通しが得られました。

またこの方程式の解はポテンシャル項
$V(r) = \int\frac{e^2}{4\pi\epsilon_0|\boldsymbol{r}-\boldsymbol{r}'|}|\varphi(|\boldsymbol{r}'|)|^2 d\boldsymbol{r}'$
が加わっているため水素原子の場合とは異なる関数形になりますが、一体のポテンシャル部分は「核からの引力クーロンポテンシャル＋雲のように広がった斥力クーロンポテンシャル」という形であり、「電子雲からの斥力で弱まった引力クーロンポテンシャル」のようなモノになりそうです。

したがって、解となる波動関数も大体が水素原子の場合に求めた固有状態、「1s状態」とか「3d状態」のような状態に近くなりそうです。

最後に、少し違う方法で上記の予想を確かめてみてHe原子の項を終わりにします。

### 遮蔽されたクーロンポテンシャルとみなす近似

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
R_{10}(r) = \left( \frac{Z}{a_0}  \right)^{3/2}2e^{-Zr/a_0} \\
Y_0^0 = \frac{1}{\sqrt{4\pi}} 
$$

を使います。


「原子核からのクーロンポテンシャルが弱まった状態」を予想しているので、原子核の電荷に対応する$Z$を変分法で最適化するパラメータと考え、変分計算を具体的に行っていきます。

ここで注意点としては、ハミルトニアンはHe原子のハミルトニアン、つまり$Z=2$として考える必要があります。

$$
\mathcal{H} = \frac{-\hbar^2}{2m}\nabla_1{}^2 + \frac{-\hbar^2}{2m}\nabla_2{}^2 
+ \frac{-2e^2}{4\pi\varepsilon_0r_1}
+ \frac{-2e^2}{4\pi\varepsilon_0r_2}
+ \frac{e^2}{4\pi\varepsilon_0r_{12}}
$$


1s状態は原点からの距離$|\boldsymbol{r}|=r$のみの関数なので以降$\varphi_{1s}(r)$と書きます。

エネルギー期待値$\langle E \rangle$は、スピン部分を先ほどと同じく先に積分して

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

です。

これを最小化するパラメータ$Z$を考えるために、$Z$の関数で書くことを目指して積分していきます。


初めに簡単な一体部分から考えていきます。^[水素様原子の一体のハミルトニアンに対して水素様原子の固有関数を考えているので、一見単純に固有値を使えばいいように思えますが（というか私が最初そのように勘違いして計算して答えが合わなかったのですが）、ハミルトニアンは$Z=2$、試行関数は$Z$を変数として扱っているので実は固有関数にならず、ちゃんと計算しないといけないのでした。]



まず運動エネルギー部分は、極座標表示のラプラシアンの$e^{-Zr/a_0}$への作用

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

この積分は以下の関係

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

を考えると実行できて、


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
    \frac{\hbar^2}{m}\left(\frac{Z}{a_0}\right)^2
\end{align*}
$$

ここに

$$
a_0 = \frac{4\pi\varepsilon_0\hbar^2}{me^2}
$$

を代入して、最終的に運動エネルギー部分は

$$
2\int\varphi_{1s}^*(r)\left( 
    \frac{-\hbar^2}{2m}\nabla{}^2 
    \right)
    \varphi_{1s}(r)d\boldsymbol{r}
=
\frac{\hbar^2}{m}\left(\frac{Z}{a_0}\right)^2
=
\frac{e^2}{4\pi\varepsilon_0}\frac{Z^2}{a_0}
$$

となります。


続いてポテンシャル部分について計算していきます。
角度部分を先に積分して、

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

先ほどと同様に定積分を計算すると、

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

というわけで先ほどの結果と合わせて、一体部分は

$$
\begin{align*}
2\int\varphi_{1s}^*(r)\left( \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-2e^2}{4\pi\varepsilon_0r} \right)
    \varphi_{1s}(r)d\boldsymbol{r}
&=
    \frac{e^2}{4\pi\varepsilon_0}\frac{Z^2}{a_0}
    -
    \frac{4e^2}{4\pi\varepsilon_0}
    \frac{Z}{a_0}
\end{align*}
$$

となります。


続いて$\langle E \rangle$の第2項

$$
\iint \varphi_{1s}^*(r)\varphi_{1s}^*(r')\frac{e^2}{4\pi\varepsilon_0\left|\boldsymbol{r}-\boldsymbol{r}'\right|}
    \varphi_{1s}(r)\varphi_{1s}(r')d\boldsymbol{r}d\boldsymbol{r}'
$$
について考えていきます。

まずは$\boldsymbol{r}'$についての積分を行います。これは有効ポテンシャル$V(r)$に対応します^[この部分の計算は[こちらのノート](http://www.th.phys.titech.ac.jp/~muto/lectures/QMII10/QMII10_chap17.pdf)を参考にしました。感謝です。]。

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
2^2\left(\frac{Z}{a_0}\right)^3
\left\{
\int_0^r
\frac{1}{r}e^{-2Zr'/a_0}r'^2dr'

+
\int_r^{\infty}
\frac{1}{r'}e^{-2Zr'/a_0}r'^2dr'
\right\}\\

&=
2^2\left(\frac{Z}{a_0}\right)^3
\left\{
\frac{1}{r}\int_0^r
e^{-2Zr'/a_0}r'^2dr'

+
\int_r^{\infty}
e^{-2Zr'/a_0}r'dr'
\right\}
\end{align*}
$$

この積分も先ほどと同様に
定積分の計算をして

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
\begin{align*}
2\int\varphi_{1s}^*(r)\left( \frac{-\hbar^2}{2m}\nabla{}^2 +  \frac{-2e^2}{4\pi\varepsilon_0r} \right)
    \varphi_{1s}(r)d\boldsymbol{r}
&=
    \frac{e^2}{4\pi\varepsilon_0}\frac{Z^2}{a_0}
    -
    \frac{4e^2}{4\pi\varepsilon_0}
    \frac{Z}{a_0}
\end{align*}
$$

と合わせて、

$$
\begin{align*}
\langle E \rangle
&=
\frac{e^2}{4\pi\varepsilon_0}\frac{Z^2}{a_0}
    -
    \frac{4e^2}{4\pi\varepsilon_0}
    \frac{Z}{a_0}
+
\frac{e^2}{4\pi\varepsilon_0}\frac{5}{8}\frac{Z}{a_0}\\

&=
-\frac{e^2}{4\pi\varepsilon_0 a_0}\left(
    Z^2 - 4Z + \frac{5}{8}Z
    \right)
\end{align*}
$$

となりました。この期待値を最小にする$Z$は、

$$
\frac{d\langle E \rangle}{dZ} = -\frac{e^2}{4\pi\varepsilon_0 a_0}\left(
    2Z - \frac{27}{8}
    \right)
= 0
$$

より、

$$
Z = \frac{27}{16}=1.6875
$$

が得られました。
手元の[猪木川合](https://bookclub.kodansha.co.jp/product?item=0000147776)をカンニングすると、この値を使って計算した$\langle E \rangle$は実測値の約98%の精度となっており、かなり良い近似であると考えられます。

これは、$Z=1.6875$とした$1s$型の波動関数つまり「原子核の電荷$Z=2$が、もう一方の電子雲によってちょっと遮蔽されて$Z=1.687$くらいに感じている波動関数」


$$
\varphi(r) = \frac{1}{\sqrt{4\pi}}\left( \frac{Z}{a_0}  \right)^{3/2}2e^{-Zr/a_0} 
$$

がHe原子中の電子の良い近似的な波動関数になっているということです。


## ここまでのまとめ
というわけでいくつかの面から、多電子原子中の電子の状態について具体的に見てきました。その結果

- 多電子原子中の電子間の相互作用は、個別の電子に着目したうえで、それ以外の電子から受ける相互作用を「電荷密度が存在確率に応じて広がったポテンシャル」（平均場ポテンシャル）として一体の相互作用として取り入れることで近似的に一体の電子の問題として扱うことができる
- そのように考えた場合、近似的に「原子核からの（引力）クーロンポテンシャルが、電子間の平均場によってちょっと遮蔽された（小さくなった）状態として有効ポテンシャルを考えることができる
- その結果、多電子原子中の電子の固有状態も、水素様原子中の電子状態と大体同じような関数で考えることができる

という雰囲気がつかめたかと思います。

2電子の問題を扱った本節に比べて一般の多電子の場合を考えるともう少し複雑にはなるのですが、大体同じような考え方・結論となります。



# 一般の多電子原子

ここまでHe原子の例を使って、かなり具体的な計算をしてきました。一般の多電子原子中の電子の有効一体ポテンシャルの考え方も、おおむね同じように、「ほかの電子から感じる斥力ポテンシャルを上手いこと平均したもの」という形で求めることができます。

最後に簡単に一般の場合を考えて本章を終わりにします。

ただ、多電子原子中の電子が感じる有効ポテンシャル（$\simeq$固体中の電子が感じる有効ポテンシャル）の近似方法については長く研究されており、現在研究で使う場合にはもっと色々と工夫がなされているようです。そこに突っ込んでいくと到底終わらなさそうなので本稿ではHF近似にとどめておくことにします。

## Hartree-Fock近似

一般の多電子原子において一つの電子が感じる有効1体ポテンシャルを考える場合も、以下のように、He原子で行った場合と同じように変分法を使います。
このような近似方法をHartree-Fock近似（HF近似と書いたりします）と呼びます。

すなわち、多変数関数である$2N$電子の波動関数$\Phi(\tau_1, \tau_2,\cdots,\tau_{2N})$^[スピン上向きと下向きで一つずつ状態を使うことから、便宜上$2N$電子とします]
を、一つのスレーター行列式$\Phi_{\rm{HF}}$を試行関数として考え、近似的に

$$
\begin{align*}
\Phi(\tau_1, \tau_2,\cdots,\tau_{2N}) &\simeq
\Phi_{\rm{HF}}(\tau_1, \tau_2,\cdots,\tau_{2N}) \\
&= 
\frac{1}{\sqrt{2N!}}
\begin{vmatrix}
\varphi_\lambda(\tau_1) & \varphi_\mu(\tau_1) & \cdots & \varphi_\xi(\tau_1)\\
\varphi_\lambda(\tau_2) & \varphi_\mu(\tau_2) & \cdots & \varphi_\xi(\tau_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\tau_{2N}) & \varphi_\mu(\tau_{2N}) & \cdots & \varphi_\xi(\tau_{2N})
\end{vmatrix}\\
&=\frac{1}{\sqrt{2N!}}\sum_P(-1)^P\hat{P}\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{2N})\\
&\equiv
 \left| \varphi_a \overline{\varphi}_a\cdots \varphi_n \overline{\varphi}_n\right|
\end{align*}
$$

と仮定します^[ちなみに、すべての組み合わせを考慮したスレーター行列式の線形結合を考えると、任意の多電子原子の波動関数を表すことができます。HF近似はこの線形結合を1つのスレーター行列式のみで近似する考え方といえます。]。[スピンを考慮した多電子原子](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)の章で断ったように、スピン状態も含んだ状態をギリシャ文字でラベルし、軌道関数のみの状態をアルファベットでラベルしています。また、軌道関数のラベルについて$a$から始めて順番に$n$で終わると$N$個にはならないのですが、そこは間に色々入れて最後のアルファベットに$n$を持ってきているという雰囲気で上手く解釈してください。

また、スレーター行列式に含まれる一体の固有関数はすべて規格直交しているものと仮定し、
このように仮定した試行関数を用いてエネルギー期待値を最小にする条件を求めます。

多電子原子のハミルトニアンは、He原子の場合のN電子版を考えて

$$

\mathcal{H} = \sum_i
\left\{\frac{-\hbar^2}{2m}\nabla_i{}^2 
+ \frac{-Ze^2}{4\pi\varepsilon_0r_i}
\right\}
+ 
\frac{1}{2}\sum_i\sum_{j\neq i} 
\frac{e^2}{4\pi\varepsilon_0|r_{i}-r_j|}\\

\equiv

\sum_i
\hat{H}_i
+ 
\frac{1}{2}\sum_i\sum_{j\neq i} 
V(\boldsymbol{r}_i,\boldsymbol{r}_j)
$$

となります。ここで$\boldsymbol{r}_j$に作用する一体のハミルトニアンを$\hat{H}_i$、二体のポテンシャルを$V(\boldsymbol{r}_i,\boldsymbol{r}_j)$と置きました。また、最後の$i,j$に関する和は二重に数えた分を$2$で割っています。

ハミルトニアンの期待値はスレーター行列式の積分

$$
\begin{align*}
&\int \Phi_{\rm{HF}}^*\mathcal{H}\Phi_{\rm{HF}}d\tau_1\cdots d\tau_{2N}

\\
&\int
\frac{1}{\sqrt{2N!}}\sum_P(-1)^P\hat{P}\varphi^*_\lambda(\tau_1)\varphi^*_\mu(\tau_2)\cdots\varphi^*_\xi(\tau_{2N})\\


&\>\>\>\>\>\>\mathcal{H}
\frac{1}{\sqrt{N!}}\sum_{P'}(-1)^{P'}\hat{P}'\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{2N})
d\tau_1\cdots d\tau_{2N}
\end{align*}
$$

となります。ここで、同じ形が$2N!$個出てくるので左側の置換を取り外して

$$
\begin{align*}
\int
\varphi^*_\lambda(\tau_1)\varphi^*_\mu(\tau_2)\cdots\varphi^*_\xi(\tau_{2N})
\mathcal{H}
\sum_{P'}(-1)^{P'}\hat{P}'\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{2N})
d\tau_1\cdots d\tau_{2N}
\end{align*}
$$

とできます。^[いつか証明を追記します]

$\mathcal{H}$を一体部分と二体部分に分けて考えていきます。

### 一体部分の積分

一体部分の積分は以下のようになります。

$$
\begin{align*}
\int
\varphi^*_\lambda(\tau_1)\varphi^*_\mu(\tau_2)\cdots\varphi^*_\xi(\tau_{2N})
\sum_i\hat{H}_i
\sum_{P'}(-1)^{P'}\hat{P}'\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{2N})
d\tau_1\cdots d\tau_{2N}\\

=
\sum_i
\int
\varphi^*_\lambda(\tau_1)\varphi^*_\mu(\tau_2)\cdots\varphi^*_\xi(\tau_{2N})
\hat{H}_i
\sum_{P'}(-1)^{P'}\hat{P}'\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{2N})
d\tau_1\cdots d\tau_{2N}

\end{align*}
$$

ここで、左側の変数$\tau_i$の関数$\varphi_\eta(\tau_i)$と、右側の置換後の変数$\tau_i$の関数$\varphi_\delta(\tau_i)$が異なる場合、異なる関数が直交しているという仮定から積分はゼロになります。従って無置換の積分だけが生き残り、「一体のハミルトニアンの各一体固有関数の期待値の総和」

$$
\sum_\nu\int \varphi_\nu^*(\tau)\hat{H}\varphi_\nu(\tau)d\tau
\equiv
\sum_\nu \langle \nu\left|\hat{H}\right|\nu\rangle
$$

となります。最後の$\langle \nu\left|\hat{H}\right|\nu\rangle$はブラケット記法と呼ばれる書き方ですが、さしあたりブラケット記法については深く触れません。「状態$\varphi^*_\nu$と$\varphi_\nu$で$\hat{H}$を挟んだ積分」くらいの意味と理解していただければOKです。

### 二体部分の積分

二体部分の積分は以下のようになります。

$$
\begin{align*}
\int
\varphi^*_\lambda(\tau_1)\varphi^*_\mu(\tau_2)\cdots\varphi^*_\xi(\tau_{2N})

\frac{1}{2}
\sum_i\sum_{j\neq i}

V(\boldsymbol{r}_i,\boldsymbol{r}_j)
\sum_{P'}(-1)^{P'}\hat{P}'\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{2N})
d\tau_1\cdots d\tau_{2N}\\

=
\frac{1}{2}
\sum_i\sum_{j\neq i}
\int
\varphi^*_\lambda(\tau_1)\varphi^*_\mu(\tau_2)\cdots\varphi^*_\xi(\tau_{2N})

V(\boldsymbol{r}_i,\boldsymbol{r}_j)

\sum_{P'}(-1)^{P'}\hat{P}'\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{2N})
d\tau_1\cdots d\tau_{2N}

\end{align*}
$$

二体部分でも、$\tau_i,\tau_j$に関する積分とそれ以外に分けて考えればよいのですが、先ほどと比べて少し複雑になります。

左側の$\tau_i, \tau_j$を変数として持つ固有関数をそれぞれ、$\varphi^*_\nu(\tau_i), \varphi^*_\eta(\tau_j)$とすると、
右側の関数が無置換の場合の$\varphi_\nu(\tau_i), \varphi_\eta(\tau_j)$か、$i$と$j$が置換された$\varphi_\nu(\tau_j), \varphi_\eta(\tau_i)$である場合以外は直交した関数の積分が出てきてしまい積分がゼロになります。


一回置換すると$(-1)$がかかることから、二体部分の積分は

$$
\begin{align*}
&\frac{1}{2}\sum_{\nu}\sum_\eta
\left\{
\iint\varphi^*_\nu(\tau) \varphi^*_\eta(\tau')
V(\boldsymbol{r}_i,\boldsymbol{r}_j)
\varphi_\nu(\tau) \varphi_\eta(\tau')
d\tau d\tau'

-
\iint\varphi^*_\nu(\tau) \varphi^*_\eta(\tau')
V(\boldsymbol{r}_i,\boldsymbol{r}_j)
\varphi_\eta(\tau) \varphi_\nu(\tau')
d\tau d\tau'

\right\}\\

&\equiv
\frac{1}{2}\sum_{\nu}\sum_\eta
\left\{
    \langle 
        \nu \eta
        \left| V(\boldsymbol{r},\boldsymbol{r}') \right|
        \nu \eta
    \rangle

    -
    \langle 
        \nu \eta
        \left| V(\boldsymbol{r},\boldsymbol{r}') \right|
        \eta \nu
    \rangle

\right\}

\end{align*}
$$

となります。最後のブラケットは、先ほどと同じような略記法です。第一項と第二項で右側の順番が異なっていることに注意してください。また和については、$\nu = \eta$の場合は第1項と第2項が打ち消しあって自動的にゼロとなるので、あえて$\nu \neq \eta$の条件はつけていません。

ここで、状態$\nu$と状態$\eta$のスピン状態のみを取り出して考えてみると、同じスピン状態を持つ場合のみ積分の第2項が値を持ち、スピン状態が異なる場合は第2項は直交したスピン関数の積分が現れ、ゼロになります。
He原子の場合は、初めからスピン状態が異なる状態をそれぞれ一体の固有関数として仮定していたのでこの項が出てこなかったわけでした。

このようにして得られたハミルトニアンの期待値に対して、He原子の場合と同じように変分法を考えることで^[本当はもう少し工夫が必要]
軌道部分の1電子波動関数$\varphi$が従う方程式

$$
\left[
\hat{H} + 2\sum_{l=a}^n \int \varphi_l^*(\boldsymbol{r}')V(\boldsymbol{r},\boldsymbol{r}')\varphi_l(\boldsymbol{r}')
d\boldsymbol{r}'
\right] \varphi(\boldsymbol{r})\\
-
\sum_{l=a}^n\left[
    \int \varphi_l^*(\boldsymbol{r}')V(\boldsymbol{r},\boldsymbol{r}')\varphi(\boldsymbol{r}')
    d\boldsymbol{r}'
\right] \varphi_l(\boldsymbol{r}) = \varepsilon \varphi(\boldsymbol{r})
$$

を得ることができます。これをHartree-Fock方程式やFock方程式などと呼びます。左辺第1項はHe原子で考えた方程式と似た形をしています。このうち

$$
\int \varphi_l^*(\boldsymbol{r}')V(\boldsymbol{r},\boldsymbol{r}')\varphi_l(\boldsymbol{r}')
    d\boldsymbol{r}'

=
\int V(\boldsymbol{r},\boldsymbol{r}')\left|\varphi_l(\boldsymbol{r}')\right|^2
    d\boldsymbol{r}'
$$

は、$l$番目の電子が作る平均場（電荷雲）ポテンシャルで、それを$l = a$から$n$まで足し合わせて2倍しているのは、スピンが異なる電子がそれぞれ2つずつあるので、全部の$2N$電子が作る平均場ポテンシャルに対応します。

この項をCoulomb演算子（クーロン演算子）と呼び、またこの形の積分をCoulomb積分と呼びます。

一方

$$
\sum_{l=a}^n\left[
    \int \varphi_l^*(\boldsymbol{r}')V(\boldsymbol{r},\boldsymbol{r}')\varphi(\boldsymbol{r}')
    d\boldsymbol{r}'
\right] \varphi_l(\boldsymbol{r}) 
$$

の部分は少し変わっていて、求めるべき関数$\varphi$が積分の中に入ってしまっています。これは積分演算子と呼ばれるそうです。
また、先ほど見たようにこの項はスピンが等しい電子間のみに働くので、$l=a$から$n$のみの$N$項のみ存在します。

この項を交換演算子、この形の積分を交換積分、交換相互作用などと呼びます。

この方程式もHe原子の場合と同様に、求めたい関数が演算子の中に入っているSelf-consistent方程式の形になっており、先述のように繰り返し計算により数値的に解くこととなります。

また、交換積分の項を積分演算子からポテンシャルの形へ近似する方法が考えられており、そのように近似した方法をHartree-Fock-Slater方程式などと呼ぶようです。
この辺りは深く突っ込むことは避けます。例えば[こちらのページ](http://cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/atom.html)に解説があります。

### ポテンシャルの球対称化近似

詳細は置いておいて、交換演算子部分を積分演算子から平均場ポテンシャルの形に近似できたことにすると、一般の多電子原子中の1電子が従う方程式は以下のような形になります。

$$
\left[
\hat{H} +2\sum_l V_l^c(\boldsymbol{r})-\sum_l V_l^{ex}(\boldsymbol{r})

\right]\varphi(\boldsymbol{r}) = \varepsilon \varphi(\boldsymbol{r})
$$



この際、He原子の場合と異なり、一般の一体固有関数$\varphi_l$は球対称ではないのでそれらが作る平均場ポテンシャルも球対称にはなっておりません。
ここからさらに、$V(\boldsymbol{r})$を角度部分で平均化することで近似的に球対称ポテンシャルと考えれば、水素原子の章で行ったように動径部分の微分方程式と角度部分の微分方程式に分解することができるようになります。

具体的には、位置$\boldsymbol{r}$のポテンシャル$V(\boldsymbol{r})$について、角度部分について積分した後に表面積$4\pi r^2$で割ることで、

$$
V(r) = \frac{1}{4\pi r^2}\int V(\boldsymbol{r})r^2\sin\theta d\theta d\phi\\
= \frac{1}{4\pi}\int V(\boldsymbol{r})\sin\theta d\theta d\phi
$$

一体の電子が感じる有効ポテンシャルを球対称な形に近似します。

このように置くことで水素様原子のシュレディンガー方程式のポテンシャル部分に、球対称な平均場ポテンシャルが付け加わった形

$$
\left(-\frac{\hbar^2}{2m} \nabla^2  -\frac{Ze^2}{4\pi\epsilon _0r}  
 + \sum_lV_l^c(r) + \sum_lV_l^{ex}(r)

\right) 

\varphi(\boldsymbol{r}) =\epsilon \varphi (\boldsymbol{r} )
$$

となり、水素原子の場合で行ったように極座標表示に変換して変数分離をすることで、3つの微分方程式

$$
-\frac{\hbar^2}{2m} \left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        -\frac{\lambda }{r^2} R
        \right)
        +
        \left(
        -\frac{Ze^2}{4\pi\epsilon _0r}
        + \sum_lV_l^c(r) + \sum_lV_l^{ex}(r)
        \right)R = \epsilon R
       ,\\
\frac{1}{\sin\theta } \frac{d}{d\theta } \left( \sin\theta \frac{d\Theta }{d\theta }   \right) + \left( \lambda - \frac{m^2}{\sin^2\theta }  \right) \Theta =0,\\
\frac{d^2\Phi }{d\phi ^2} + m^2\Phi =0.
$$

へ帰着させることができます。

角度部分は水素様原子とまったく同じ形なので、固有関数は球面調和関数

$$
Y_l^m(\theta,\phi)
$$

で書くことができ、動径関数のみを「大体クーロンポテンシャル」の微分方程式をself-consistentに解くことで近似的な一体の固有関数を得られることがわかりました。すなわち、多電子原子中の電子であっても、近似的に水素様原子の場合と同様$1s, 2p, 3d,\cdots$状態のようなものだと考えられるということです。
今後、特にtight-binding近似を考える際に多電子原子の電子の状態を「$s$状態」や「$p$状態」などと表しますが、上記のような、大体水素様原子の固有状態で表せる、という仮定のもと使用しているものと考えてください。

# まとめ

He原子のまとめの項のほぼ再掲になりますが、本章を通して多電子原子中の電子状態や、電子が感じる有効ポテンシャルについて以下のことがわかりました。

- 多電子原子中の電子間の相互作用は、個別の電子に着目したうえで、それ以外の電子から受ける相互作用を「電荷密度が存在確率に応じて広がったポテンシャル」（平均場ポテンシャル）として一体の相互作用として取り入れることで近似的に一体の電子の問題として扱うことができる
- そのように考えた場合、近似的に「原子核からの（引力）クーロンポテンシャルが、電子間の平均場によってちょっと遮蔽された（小さくなった）状態として有効ポテンシャルを考えることができる
- その結果、多電子原子中の電子の固有状態も、水素様原子中の電子状態と大体同じような関数で考えることができる

現在はHF近似よりもさらに進んだ理論が色々とあり、本章に書いた内容だけでは不十分ではあるのですが、おおむね有効ポテンシャルの考え方について理解できたことにしておき、いずれ余裕ができたらそのあたりも学んでいければと思います。

この内容を踏まえ、次章からいよいよ原子核を増やした場合の性質、つまり固体物理の内容に入っていきます。その際はほとんどの場合、固体を構成する原子の価電子に着目し、内殻電子との相互作用を上記のような一体有効ポテンシャルとして近似する扱いをします。
以降、価電子が感じる有効ポテンシャルの関数形は明記しませんが、本章で扱ったような考えのもと、数値的に得られていると考えて進めていきます。
