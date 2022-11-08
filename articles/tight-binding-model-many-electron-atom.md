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

## Self-consistent近似

上記電子間相互作用を無視した解を求めましたが、実際には電子はそれぞれもう一方の電子からの斥力ポテンシャルを感じているはずです。その効果を取り入れるために、
参考：http://cms.phys.s.u-tokyo.ac.jp/~naoki/CIPINTRO/CIP/atom.html

# 一般の多電子原子
## Hartree近似

## Hartree-Fock近似

## 多電子原子内の電子の波動関数と平均場ポテンシャル

