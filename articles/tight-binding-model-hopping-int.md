---
title: "飛び移り積分（Hopping Integral）の物理的な意味・Wannier関数の従う方程式"
emoji: "🐸"
type: "idea" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: false
---
# はじめに

前章でTight-bindingモデルを導出していく中で、原子軌道関数や局所ポテンシャルの積分で定義された重なり積分

$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r} \simeq \delta_{m,n}\delta_{\boldsymbol{R},\boldsymbol{R}'}
$$

や飛び移り積分

$$
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_I)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d\boldsymbol{r}
 \equiv
 - t_{\boldsymbol{R}_I}^{n,m}
$$

が定義され、それぞれ重なり積分は「小さな値」として近似されたり、飛び移り積分はバンドの形状に大きな影響を与えたり、Tight-bindingモデルの構築に当たり重要な役割を果たしていました。

それぞれの積分は、「何を計算すればよいか」は式を見れば明らかなのですが、じゃあ実際どのような物理的なイメージと対応しているのかと考え始めてみると分かりそうでよくわかりません（でした、私は）。
ググってみてもイマイチ自分が欲しい答えも見つかりません（ググり力が足りないだけかもしれません）。

というわけで、色々考えて自分なりに納得した内容について、本章で飛び移り積分について、次の章重なり積分について、書き残しておこうと思います。
**間違った内容を書いてたりする可能性大なので注意してください。**

# 孤立原子軌道に対する飛び移り積分について

まずは比較的（多分）まともなことが書けそうな飛び移り積分の方から始めていきます。

「飛び移り積分」と名前がついている通り、格子点間の電子の移動に関係が深い積分であることが予測され、さらに少し先取りして書くと、第二量子化表示をした場合


$$
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }-t_{ij}a^\dagger_{i\sigma }a_{j\sigma } + h.c.
$$

のように書かれ、「サイトjの電子が、サイトiの電子に変わる」＝飛び移っていくようなそんな式の係数になっていたりします。

というわけで何か電子の移動確率みたいなもんなのかな？という印象を強く受けるわけですが、具体的にどのような物理量や物理的現象に対応しているのか、分かるようでこれを見ただけではわかりません。「あるサイトの原子軌道と、隣のサイトの局所ポテンシャルと原子軌道との積を積分することはわかったけどそれでその値は、何を意味するの？」って感じです。

というわけで順を追って、飛び移り積分と物理的な現象？を結びつけていきたいと思います

## 原子軌道関数の従う方程式

まず初めに、前章と同じく固体のポテンシャルを1電子近似し、かつ局所ポテンシャルの和で近似した1電子ハミルトニアン

$$
\hat{H}^{\rm c} = \frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
$$

を考えます（ハミルトニアンの添え字"c"はCrystalの"c"のつもり）。ここで、固体の3方向の基本併進ベクトルを$\boldsymbol{a}_i, i=1,2,3$、格子ベクトルを$\boldsymbol{R} = \sum_in_i\boldsymbol{a}_i, n_i = 0,\pm 1,\pm 2\cdots$として、$\boldsymbol{a}_i$方向に$N_i\boldsymbol{a}_i$進むと元に戻る周期的境界条件を考え、局所ポテンシャルの総和$\sum_{\boldsymbol{R}}$も周期的境界条件の1周期分（≒結晶全体）について和を取ります。

また、前章と同じように原点のポテンシャルだけ抜き出した孤立原子ハミルトニアン$\hat{H}^{\rm a} = \frac{-\hbar^2}{2m}\nabla^2 + V(\boldsymbol{r})$に対して、固有値（原子準位）$\varepsilon_m^{\rm a}$を持つ固有状態である原子軌道を$\phi_m(\boldsymbol{r})$と書きます。

ここである原子軌道$\phi_m$に、結晶のハミルトニアンをおもむろに掛けてみますと、

$$
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
\right) \phi_m(\boldsymbol{r})
$$

固有関数$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$がBloch和$\Phi_{m,\boldsymbol{k}}(\boldsymbol{r}) = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r} - \boldsymbol{R})$で展開できたことから、原子軌道$\phi_m(\boldsymbol{r})$にハミルトニアンが作用した後の関数も、関数系$\{\phi_m(\boldsymbol{r} - \boldsymbol{R})\}$で

$$
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
\right) \phi_m(\boldsymbol{r}) = \sum_{m,\boldsymbol{R}}c_{m,\boldsymbol{R}}\phi_{m}(\boldsymbol{r} - \boldsymbol{R})
$$

と展開できます。

:::details 証明っぽい内容

突然原子軌道関数を格子点だけ平行移動させた関数が登場しましたが、本当にこの関数系で任意の関数を展開できるのでしょうか？

これは以下のように考えれば良さそうです^[もっと正当な証明があるのだろうとは思いますが、あまりここにこだわってもしょうがないですし、よくわかっていないのでオレオレ証明で一旦ごまかしておきます。教科書を眺めていると「Sturm–Liouvilleの理論」のようなキーワードがチラ見えしますがまあ置いておきます。]：

まず、何度も見てきたように、結晶のハミルトニアン$\hat{H}^{\rm c}$の固有値$\varepsilon_{n,\boldsymbol{k}}$を持つ固有状態$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$：

$$
\hat{H}^{\rm c} \varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})
$$

は、孤立原子ハミルトニアン$\hat{H}^{\rm a} = \frac{-\hbar^2}{2m}\nabla^2 + V(\boldsymbol{r})$に対して、固有値（原子準位）$\varepsilon_m^{\rm a}$を持つ固有状態である原子軌道$\phi_m(\boldsymbol{r})$からなるBloch和$\Phi_{m,\boldsymbol{k}}(\boldsymbol{r}) = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r} - \boldsymbol{R})$で

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) =
\sum_m b_m^n\Phi_m(\boldsymbol{r}) ,\\

\Phi_{m,\boldsymbol{k}}(\boldsymbol{r}) = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{m}(\boldsymbol{r} - \boldsymbol{R})
$$

と展開できます。行列とベクトルっぽい物で表現すると、

$$
(B)_{nm} = b_{m}^{n}
$$

から行列$B$を定義して、

$$

\begin{bmatrix}
\varphi_{1,\boldsymbol{k}}(\boldsymbol{r}) \\
\varphi_{2,\boldsymbol{k}}(\boldsymbol{r}) \\
\vdots\\
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) \\
\vdots
\end{bmatrix}

=
B
\begin{bmatrix}
\Phi_{1,\boldsymbol{k}}(\boldsymbol{r}) \\
\Phi_{2,\boldsymbol{k}}(\boldsymbol{r}) \\
\vdots\\
\Phi_{n,\boldsymbol{k}}(\boldsymbol{r}) \\
\vdots
\end{bmatrix}

$$

のような変換を考えることができます。


一方この展開を逆に考えると、$B^{-1}B = I$（$I$は単位行列）を満たす逆行列$B^{-1}$の要素を

$$
(B^{-1})_{nm} \equiv \tilde{b}_{m}^{n}
$$

と書くことにして、逆変換

$$

\begin{bmatrix}
\Phi_{1,\boldsymbol{k}}(\boldsymbol{r}) \\
\Phi_{2,\boldsymbol{k}}(\boldsymbol{r}) \\
\vdots\\
\Phi_{n,\boldsymbol{k}}(\boldsymbol{r}) \\
\vdots
\end{bmatrix}
=
B^{-1}


\begin{bmatrix}
\varphi_{1,\boldsymbol{k}}(\boldsymbol{r}) \\
\varphi_{2,\boldsymbol{k}}(\boldsymbol{r}) \\
\vdots\\
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) \\
\vdots
\end{bmatrix}

$$

を定義できます。



これらとBloch和の「逆変換」$\phi_m(\boldsymbol{r}-\boldsymbol{R}) = N^{-1}\sum_{\boldsymbol{k}}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})$を用いて、原子軌道関数を固有関数の線形結合で表す以下の関係式が得られます：

$$
\phi_{m}(\boldsymbol{r}-\boldsymbol{R}) = \frac{1}{N}\sum_{n,\boldsymbol{k}}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\tilde{b}_m^n\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}).
$$

特に、$\boldsymbol{R} = \boldsymbol{0}$の関係式は以下の通りです：

$$
 \phi_{m}(\boldsymbol{r}) = \frac{1}{N}\sum_{n,\boldsymbol{k}}\tilde{b}_m^n\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}).
$$

さて、この関係式を初めの式へ代入すると、$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$は$\hat{H}^{\rm c}$の固有関数なので、


$$
\begin{align*}
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
\right) \phi_m(\boldsymbol{r}) 

&=
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
\right) 
\frac{1}{N}\sum_{n,\boldsymbol{k}}\tilde{b}_m^n\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})
\\

&=
\frac{1}{N}\sum_{n,\boldsymbol{k}}\varepsilon_{n,\boldsymbol{k}}\tilde{b}_m^n\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})

\end{align*}
$$

となりますが、さらに今度は$\varphi_{n\boldsymbol{k}}$を原子軌道関数（Bloch和）で展開して、

$$
\begin{align*}
（上式右辺）&=
\frac{1}{N}\sum_{n,\boldsymbol{k}}\varepsilon_{n,\boldsymbol{k}}\tilde{b}_m^n\sum_{m'} b_{m'}^n\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{m'}(\boldsymbol{r} - \boldsymbol{R})
\end{align*}
$$

さらにBloch波数の関数$\varepsilon_{n,\boldsymbol{k}}$のFourier展開（Bloch関数に加え、固有関数もBloch波数について逆格子ベクトルを周期とする周期関数なのでした。というわけで、Bloch関数と同様に、直接格子の格子ベクトルを用いてFourier展開できます）

$$
\varepsilon_{n,\boldsymbol{k}} = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{R}},\\
\varepsilon_{n,\boldsymbol{R}} = \frac{1}{N}\sum_{\boldsymbol{k}}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{k}} 
$$

を代入して、

$$
\begin{align*}
（上式右辺）&=
\frac{1}{N}\sum_{n,\boldsymbol{k}}\sum_{\boldsymbol{R}'}e^{i\boldsymbol{k}\cdot\boldsymbol{R'}}\varepsilon_{n,\boldsymbol{R}'}\tilde{b}_m^n\sum_{m'} b_{m'}^n\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_{m'}(\boldsymbol{r} - \boldsymbol{R})



\end{align*}
$$

となります。ここで$\boldsymbol{k}$の総和が計算できるので、$\sum_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R} + \boldsymbol{R}')} = N\delta_{\boldsymbol{R},-\boldsymbol{R}'}$を用いて

$$
\begin{align*}
（上式右辺）&=
\sum_{m'}
\sum_{\boldsymbol{R}}

\left(\sum_{n}\varepsilon_{n,-\boldsymbol{R}}\tilde{b}_m^n b_{m'}^n
\right)\phi_{m'}(\boldsymbol{r} - \boldsymbol{R})


\\
&\equiv
\sum_{m',\boldsymbol{R}}c_{m',\boldsymbol{R}}^{m}\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}).
\end{align*}
$$

以上のように、係数の具体的な中身はともかく、「原子軌道関数に結晶のハミルトニアンを作用させた後の関数」が、関数系$\{\phi_m(\boldsymbol{r}-\boldsymbol{R})\}$で展開できることを示せました。

:::

ここで、前章に引き続きTight-bindingな近似を取り入れ、異なる原子軌道間の重なり積分をゼロと置く：


$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r} \simeq \delta_{m,n}\delta_{\boldsymbol{R},\boldsymbol{R}'}
$$

と、展開係数$c_{m',\boldsymbol{R}}^m$は、$\phi_{m'}(\boldsymbol{r} - \boldsymbol{R})$の係数が$c_{m',-\boldsymbol{R}}^m$であることに注意して、

$$
c_{m',\boldsymbol{R}}^m = \int\phi^*_{m'}(\boldsymbol{r} + \boldsymbol{R})\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}'}V(\boldsymbol{r}-\boldsymbol{R}')
\right)\phi_m(\boldsymbol{r})d\boldsymbol{r}
$$

となります。さらに、孤立原子ハミルトニアン部分の固有関数であることを利用して、

$$
\begin{align*}
（上式右辺）&=

\varepsilon_m^{\rm a}\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R})\phi_m(\boldsymbol{r})d \boldsymbol{r} + 
\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R})
\left(
  \sum_{\boldsymbol{R}' \neq -\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}')
\right)
\phi_m(\boldsymbol{r})d \boldsymbol{r} \\

&=
\varepsilon_m^{\rm a}\delta_{m,m'}\delta_{\boldsymbol{R},\boldsymbol{0}}+ 
\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R})
\left(
  \sum_{\boldsymbol{R}' \neq -\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}')
\right)
\phi_m(\boldsymbol{r})d \boldsymbol{r} 

\end{align*}
$$

と書き換えることができます。ついに「飛び移り積分」っぽい形の積分がでてきましたね。

ここで$\boldsymbol{R} = \boldsymbol{0}$とその他で分けると、

### $\boldsymbol{R} = \boldsymbol{0}$

$$
\begin{align*}
c_{m',\boldsymbol{0}}^m &\simeq

\varepsilon_m^{\rm a}\delta_{mm'}

+
\int\phi_{m'}^*(\boldsymbol{r})
\left(
  \sum_{\boldsymbol{R}' \neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}')
\right)
\phi_m(\boldsymbol{r})d \boldsymbol{r} \\

&=
\varepsilon_m^{\rm a}\delta_{mm'}

+\Delta\varepsilon_{m'm}
\end{align*}
$$

と、原子準位$\varepsilon_m^{\rm a}$と結晶場積分$\Delta\varepsilon_{m'm}$で書けることがわかります。

### $\boldsymbol{R} \neq \boldsymbol{0}$

また$\boldsymbol{R} \neq \boldsymbol{0}$の場合は、前章と同様にTight-binding近似として、（最）隣接格子間の飛び移り積分ではないものをゼロと置くと、


$$
c_{m',\boldsymbol{R}}^m \simeq

\delta_{\boldsymbol{R},\boldsymbol{R}_I}\sum_{\boldsymbol{R}_I} 
\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R}_I)
  V(\boldsymbol{r})
\phi_m(\boldsymbol{r})d \boldsymbol{r} 
$$

を得ます。ここで積分$\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R}_I)V(\boldsymbol{r})\phi_m(\boldsymbol{r})d \boldsymbol{r}$は、ちゃんと書けば積分範囲が結晶全体の定積分ですが、積分変数を$\boldsymbol{r} +\boldsymbol{R}_I\rightarrow \boldsymbol{r}$と置きなおすと、周期的境界条件より積分範囲は変わらず、飛び移り積分の形

$$
\int_V\phi_{m'}^*(\boldsymbol{r})
  V(\boldsymbol{r}-\boldsymbol{R}_I)
\phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d \boldsymbol{r} \equiv - t_{\boldsymbol{R}_I}^{m'm}
$$

と一致します。以上より、展開係数は

$$
c_{m',\boldsymbol{R}}^m \simeq
(\varepsilon_m^{\rm a}\delta_{mm'} +\Delta\varepsilon_{m'm})\delta_{\boldsymbol{R},\boldsymbol{0}}

-
\delta_{\boldsymbol{R},\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m'm}
$$

となり、原子軌道関数に結晶のハミルトニアンを作用させた関係式は

$$
\begin{align*}
\hat{H}^{\rm c} \phi_m(\boldsymbol{r})
&=
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
\right) \phi_m(\boldsymbol{r}) \\


&\simeq
\sum_{m,\boldsymbol{R}}
\left\{
  (\varepsilon_m^{\rm a}\delta_{mm'} +\Delta\varepsilon_{m'm})\delta_{\boldsymbol{R},\boldsymbol{0}}
  -
  \delta_{\boldsymbol{R},\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{m'm}
\right\}
\phi_{m}(\boldsymbol{r} - \boldsymbol{R})\\

&=
(\varepsilon_m^{\rm a} + \Delta\varepsilon_{mm})\phi_m(\boldsymbol{r})

+
\sum_{m'} \Delta\varepsilon_{m'm}\phi_{m'}(\boldsymbol{r})\\

&\>\>\>\>+\sum_{\boldsymbol{R}_I}\sum_{m'}(-t_{\boldsymbol{R}_I}^{m'm})\phi_{m'}(\boldsymbol{r} + \boldsymbol{R}_I)

\end{align*}
$$

となることが分かりました。（$\phi_{m'}(\boldsymbol{r} -\boldsymbol{R})$の係数を$c_{m'-\boldsymbol{R}}^m$と置いたので、$\boldsymbol{R}$の符号がプラスになっています。）


・・・だからどうしたという感じですが、ここで記憶の片隅から、波動関数にハミルトニアンを作用させた意味を引っ張り出してきます。



## 波動関数の時間発展

さて、定常状態ばっかり扱っていてもはや記憶の彼方に葬り去られていましたが、本来のシュレーディンガー方程式は時間と位置に関する微分方程式

$$
i\hbar \frac{\partial}{\partial t}\psi(\boldsymbol{r},t) = \hat{H}\psi(\boldsymbol{r},t)
$$

なのでした。さらに、時刻$t$から微小時間$\Delta t$の変化を考えると、左辺の偏微分を

$$
\frac{\partial}{\partial t}\psi(\boldsymbol{r},t)\simeq

\frac{\psi(\boldsymbol{r},t + \Delta t) - \psi(\boldsymbol{r},t)}{\Delta t}
$$

として、シュレーディンガー方程式の両辺を$i\hbar$で割って整理すると、微小時間についての波動関数の時間発展の式

$$
\frac{\psi(\boldsymbol{r},t + \Delta t) - \psi(\boldsymbol{r},t)}{\Delta t}

\simeq

\frac{-i}{\hbar}\hat{H}\psi(\boldsymbol{r},t)\\

\Rightarrow
\psi(\boldsymbol{r},t + \Delta t) \simeq 

\psi(\boldsymbol{r},t ) + \frac{-i}{\hbar}\Delta t\hat{H}\psi(\boldsymbol{r},t)
$$

を得ることができます。これで準備が整ったので、**結晶の中、$t=0$で原点に一つだけ孤立原子軌道がある場合の時間発展**を考えてみましょう。すなわちハミルトニアン、$t=0$の波動関数を

$$
\begin{align*}
&\hat{H} = \hat{H}^{\rm c} 
=
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}-\boldsymbol{R})
\right),\\
&\psi(\boldsymbol{r},0) = \phi_m(\boldsymbol{r})
\end{align*}
$$

の時間発展を考えます。これは先ほど得た関係式より、


$$
\begin{align*}
\psi(\boldsymbol{r},\Delta t) &\simeq 

\psi(\boldsymbol{r},0 ) + \frac{-i}{\hbar}\Delta t\hat{H}\psi(\boldsymbol{r},0)
\\

&\simeq
\phi_m(\boldsymbol{r})
+
\frac{-i}{\hbar}\Delta t\hat{H}\phi_m(\boldsymbol{r})
\\

&\simeq

\left\{1-\frac{-i}{\hbar}\Delta t
(\varepsilon_m^{\rm a} + \Delta\varepsilon_{mm})
\right\}
\phi_m(\boldsymbol{r})\\

&\>\>\>\>+
\frac{-i}{\hbar}\Delta t\sum_{m'} \Delta\varepsilon_{m'm}\phi_{m'}(\boldsymbol{r})\\

&\>\>\>\>+\frac{-i}{\hbar}\Delta t\sum_{\boldsymbol{R}_I}\sum_{m'}(-t_{\boldsymbol{R}_I}^{m'm})\phi_{m'}(\boldsymbol{r} + \boldsymbol{R}_I)


\end{align*}
$$

が得られます。すなわち、原点にポツンと1つ原子軌道があるとき、微小時間後の波動関数を各格子点に平行移動した、様々な準位の原子軌道$\{\phi_{m'}(\boldsymbol{r} + \boldsymbol{R}\}$の重ね合わせで書くことができました。この時の$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R})$の係数はかなり雑に言うと「電子が格子点$-\boldsymbol{R}$上で、$m'$の原子軌道で存在する確率」つまり「原点の格子点$\boldsymbol{R} = \boldsymbol{0}$にいた状態$m$の電子が、格子点$-\boldsymbol{R}$に状態$m'$で「飛び移る」確率」と考えられそうです。

#### 注意点1
ここで位置を表す変数が$\boldsymbol{r}$と、$\boldsymbol{R}$の二つ出てきて紛らわしいのですが、前者$\boldsymbol{r}$は「電子の波動関数の形」$\simeq$「電子の確率雲の形」に対応する変数で、後者$\boldsymbol{R}$は「波動関数の中心」$\simeq$「固体の中の電子の場所」に対応する変数です。この表現もかなり正確性がアヤシイですが。。。

#### 注意点2

正確（？）には、波動関数は観測することができず、観測できるのはエルミート演算子で表される何かしらの物理量です。そして、ある波動関数を何らかの演算子の固有状態で展開していた時、物理量の測定を通してその固有状態である波動関数が確定します（なんかこの書き方もイマイチな気がする・・）。そこで、波動関数$\phi_m(\boldsymbol{r}-\boldsymbol{R})$は何か物理量に対応する演算子の固有状態になっているかというと、そこが良くわかっていないので、例えば上記のような時間発展を考えたとしても、時刻$\Delta t$後に、「電子が格子点$\boldsymbol{R}$上で、$m'$の原子軌道で存在する」かどうかを確認できるかというとそういうわけではないように思えます。

## 原子軌道の飛び移り積分の物理的意味


特に、微小時間$\Delta t$後に、電子が原点から（隣接）格子ベクトル$-\boldsymbol{R}_I$で指定される格子点に局在した原子軌道$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R})$となっているような確率は、飛び移り積分$t_{\boldsymbol{R}_I}^{m'm}$：

$$
\begin{align*}
t_{\boldsymbol{R}_I}^{m'm} &= 
-
\int_V\phi_{m'}^*(\boldsymbol{r})
  V(\boldsymbol{r}-\boldsymbol{R}_I)
\phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d \boldsymbol{r} \\

&=
-
\int_V\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R}_I)
  V(\boldsymbol{r})
\phi_m(\boldsymbol{r})d \boldsymbol{r} 
\end{align*}
$$

に比例することがわかります。ここで、最後の式では$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R})$と座標をそろえるために再度積分変数を変更しました。



以上のように（きわめて怪しい）議論を繰り広げることで、飛び移り積分はその名前の通り、**原点（あるいはある格子点$\boldsymbol{R}'$）にいる$m$状態の原子軌道$\phi_m(\boldsymbol{r})$（$\phi_m(\boldsymbol{r}-\boldsymbol{R}')$）が、（微小時間後に）状態$m'$になって原点（あるいはある格子点$\boldsymbol{R}'$）から格子ベクトル$-\boldsymbol{R}_I$離れた場所へ、$\phi_{m'}(\boldsymbol{r}+\boldsymbol{R}_I)$（$\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}'+\boldsymbol{R}_I)$）として飛び移って行く確率**（に比例する量）を表すことが分かりました。

ただ、飛び移り積分の（教科書通りの）定義を採用すると、「飛び移り先」が$-\boldsymbol{R}_I$になってしまうのが気になるところなのです。何か対称性をちゃんと考えれば解消されるか、どこかで勘違いしているかもしれませんが、一旦思考停止しておきます。この後来る第二量子化とかを考えていけば何か解消されるかもしれません。

~~まあ、これがわかったところで、微小時間後の波動関数が測定できるわけでもないし、何の意味があるのかは分かりませんが。というか自己満足ですが。~~






# （おまけ）Wannier関数の従うべき方程式と時間発展

最後に、~~これも何の役に立つのかわかりませんが、~~局在基底であるWannier関数の関係式と、時間発展についておまけで触れておきます。

先ほどと同様に固体のポテンシャルを1電子近似した1電子ハミルトニアン

$$
\hat{H}^{\rm c} = \frac{-\hbar^2}{2m}\nabla^2 + U(\boldsymbol{r})
$$

を考えます。ここで$U(\boldsymbol{r})$は、固体の3方向の基本併進ベクトルを$\boldsymbol{a}_i, i=1,2,3$、格子ベクトルを$\boldsymbol{R} = \sum_in_i\boldsymbol{a}_i, n_i = 0,\pm 1,\pm 2\cdots$として$U(\boldsymbol{r} + \boldsymbol{R}) =U(\boldsymbol{r})$を満たす周期ポテンシャルです。

何度も見てきたように、このハミルトニアンの固有値$\varepsilon_{n,\boldsymbol{k}}$を持つ固有状態$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$：

$$
\hat{H}^{\rm c} \varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})
$$

は、局在関数Wannier関数$w_{n,\boldsymbol{R}}(\boldsymbol{r})=w_n(\boldsymbol{r}-\boldsymbol{R})$で

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) =\sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}},\\
w_{n,\boldsymbol{R}}(\boldsymbol{r})
=
\frac{1}{N}\sum_{\boldsymbol{k}} \varphi_n(\boldsymbol{r},\boldsymbol{k})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
$$

と展開できます。このWannier関数の展開を1体シュレーディンガー方程式に代入すると、固体の1電子ハミルトニアンに対するWannier関数の満たす関係式

$$
\hat{H}^{\rm c}w_{n,\boldsymbol{R}} = \sum_{\boldsymbol{R}'}\varepsilon_{n,\boldsymbol{R}'}w_{n,\boldsymbol{R}-\boldsymbol{R}'} 
$$

を導くことができます。ここで$\varepsilon_{n,\boldsymbol{R}}$は、Bloch波数の関数$\varepsilon_{n,\boldsymbol{k}}$のFourier展開

$$
\varepsilon_{n,\boldsymbol{k}} = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{R}},\\
\varepsilon_{n,\boldsymbol{R}} = \frac{1}{N}\sum_{\boldsymbol{k}}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{k}} 
$$

です。

Wannier関数の正規直交性$\int_V w_{n'}^*(\boldsymbol{r}-\boldsymbol{R}')w_n(\boldsymbol{r} - \boldsymbol{R})dr =\delta_{n,n'}\delta_{\boldsymbol{R},\boldsymbol{R}'}$から、$\varepsilon_{n,\boldsymbol{R}} = \int w_{n,\boldsymbol{0}}^*(\boldsymbol{r})\hat{H}^{\rm c}w_{n,-\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}$です。これはBloch関数が満たすシュレーディンガー方程式に$\varepsilon_{n,\boldsymbol{k}} = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{R}}$を代入して$e^{-i\boldsymbol{k}\cdot\boldsymbol{R}'}$をかけて積分しても示すことができます

:::details 証明

まずWannier関数での展開をブロッホ関数が満たすシュレーディンガー方程式に代入する：

$$
\hat{H}^{\rm c} \sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}}= \varepsilon_{n,\boldsymbol{k}}\sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}}.
$$

次に固有値$\varepsilon_{n,\boldsymbol{k}}$について、これもBloch波数$\boldsymbol{k}$と逆格子ベクトル$\boldsymbol{K}$に関して、$\varepsilon_{n,\boldsymbol{k} + \boldsymbol{K}} = \varepsilon_{n,\boldsymbol{k}}$の周期性より、Bloch関数と同様に「逆格子ベクトルの逆格子ベクトル」格子ベクトル$\boldsymbol{R}$を用いてFourier展開できる：

$$
\varepsilon_{n,\boldsymbol{k}} = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{R}},\\
\varepsilon_{n,\boldsymbol{R}} = \frac{1}{N}\sum_{\boldsymbol{k}}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{k}} 
$$

ので、それもシュレーディンガー方程式に代入する：

$$
\hat{H}^{\rm c} \sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}}

=
 \sum_{\boldsymbol{R}'}e^{i\boldsymbol{k}\cdot\boldsymbol{R}'}\varepsilon_{n,\boldsymbol{R}'}
 \sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}}.
$$

両辺に$e^{-i\boldsymbol{k}\cdot\boldsymbol{R}''}$をかけて$\boldsymbol{k}$についてBZ内での積分（和）を取ると、$N^{-1}\sum_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R} - \boldsymbol{R}'')} = \delta_{\boldsymbol{R},\boldsymbol{R}''}$より、

$$
\begin{align*}
（左辺）&=
\hat{H}^{\rm c} \sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})N^{-1}\sum_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R} - \boldsymbol{R}'')}\\

&=
\hat{H}^{\rm c} \sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})\delta_{\boldsymbol{R},\boldsymbol{R}''}\\

&=
\hat{H}^{\rm c} w_{n,\boldsymbol{R}''}(\boldsymbol{r}).

\end{align*}
$$

$$
\begin{align*}
（右辺）&=
 \sum_{\boldsymbol{R}'} \sum_{\boldsymbol{R}}\varepsilon_{n,\boldsymbol{R}'}
w_{n,\boldsymbol{R}}(\boldsymbol{r})
N^{-1}\sum_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R} + \boldsymbol{R}' - \boldsymbol{R}'')}\\

&=
 \sum_{\boldsymbol{R}'} \sum_{\boldsymbol{R}}\varepsilon_{n,\boldsymbol{R}'}
w_{n,\boldsymbol{R}}(\boldsymbol{r})
\delta_{\boldsymbol{R},\boldsymbol{R}''-\boldsymbol{R}'}\\

&=

\sum_{\boldsymbol{R}'}\varepsilon_{n,\boldsymbol{R}'}
w_{n,\boldsymbol{R}''-\boldsymbol{R}'}(\boldsymbol{r}).

\end{align*}
$$

ここで記号を整えるために左辺、右辺ともに$\boldsymbol{R}'' = \boldsymbol{R}$と置きなおすと、等式

$$
\hat{H}^{\rm c} w_{n,\boldsymbol{R}}(\boldsymbol{r})
=
\sum_{\boldsymbol{R}'}\varepsilon_{n,\boldsymbol{R}'}
w_{n,\boldsymbol{R}-\boldsymbol{R}'}(\boldsymbol{r})
$$

が得られる。
::::


## Wannier関数の時間発展

先ほどと同様に$t=0$で原点に局在したWannier関数の時間発展を考えると、$\psi(\boldsymbol{r},0 ) =w_n(\boldsymbol{r})$として、

$$
\begin{align*}
\psi(\boldsymbol{r},\Delta t) &\simeq 

\psi(\boldsymbol{r},0 ) + \frac{-i}{\hbar}\Delta t\hat{H}\psi(\boldsymbol{r},0)
\\

&=
w_n(\boldsymbol{r})
+
\frac{-i}{\hbar}\Delta t\hat{H}w_n(\boldsymbol{r})
\\

&=

w_n(\boldsymbol{r})
+
\frac{-i}{\hbar}\Delta t\sum_{\boldsymbol{R}'}\varepsilon_{n,\boldsymbol{R}'}
w_{n,\boldsymbol{R}-\boldsymbol{R}'}(\boldsymbol{r})\\

\end{align*}
$$

となります。ここに$\Delta t$の時間発展を繰り返し考えていっても、どこまで行ってもWannier関数の関数形はそのままで、各格子点にいる確率だけが変わっていくと考えられます。

従って、結晶のポテンシャルの中で局在した一つのWannier関数は、形を変えることなく格子点を「飛び移って」行くことがわかります。

# おわりに






