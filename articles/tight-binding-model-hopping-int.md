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
**間違った内容や勘違いを書いてたりする可能性は大いにあるので注意してください。**

# 飛び移り積分について

まずは比較的（多分）まともなことがかけそうな飛び移り積分の方から始めていきます。

「飛び移り積分」と名前がついている通り、格子点間の電子の移動に関係が深い積分であることが予測され、さらに少し先取りして書くと、第二量子化表示をした場合


$$
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }-t_{ij}a^\dagger_{i\sigma }a_{j\sigma } + h.c.
$$

のように書かれ、「サイトjの電子が、サイトiの電子に変わる」＝飛び移っていくようなそんな式の係数になっていたりします。

というわけで何か電子の移動確率みたいなもんなのかな？という印象を強く受けるわけですが、具体的にどのような物理量や物理的現象に対応しているのか、分かるようでこれを見ただけではわかりません。「あるサイトの原子軌道と、隣のサイトの局所ポテンシャルと原子軌道との積を積分することはわかったけどそれでその値は、何を意味するの？」って感じです。

というわけで順を追って、飛び移り積分と物理的な現象？を結びつけていきたいと思います

## Wannier関数の従うべき方程式

まず初めに、固体のポテンシャルを1電子近似した1電子ハミルトニアン

$$
\hat{H}^{\rm c} = \frac{-\hbar^2}{2m}\nabla^2 + U(\boldsymbol{r})
$$

を考えます（ハミルトニアンの添え字"c"はCrystalの"c"のつもり）。ここで$U(\boldsymbol{r})$は、固体の3方向の基本併進ベクトルを$\boldsymbol{a}_i, i=1,2,3$、格子ベクトルを$\boldsymbol{R} = \sum_in_i\boldsymbol{a}_i, n_i = 0,\pm 1,\pm 2\cdots$として$U(\boldsymbol{r} + \boldsymbol{R}) =U(\boldsymbol{r})$を満たす周期ポテンシャルです。

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

です。^[なお、今回の話と直接関係ないですがWannier関数の正規直交性$\int_V w_{n'}^*(\boldsymbol{r}-\boldsymbol{R}')w_n(\boldsymbol{r} - \boldsymbol{R})dr =\delta_{n,n'}\delta_{\boldsymbol{R},\boldsymbol{R}'}$から、$\varepsilon_{n,\boldsymbol{R}} = \int w_{n,\boldsymbol{R}}^*(\boldsymbol{r})\hat{H}^{\rm c}w_{n,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}$です。これはBloch関数が満たすシュレーディンガー方程式に$\varepsilon_{n,\boldsymbol{k}} = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{R}}$を代入して$e^{-i\boldsymbol{k}\cdot\boldsymbol{R}'}$をかけて積分しても示すことができます]

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

ここまでは（1電子近似を除いて）近似のない一般的な関係です。（なお、この関係式は第二量子化表示を導く際にも使います。）

## 時間発展の方程式（名前？）
//こっからはちゃんと書かないといけなくなり、またそれなりのボリュームがあるけど、とりあえず↑の関係式をスッキリ証明できたので、きりがよい！

## 原子軌道での展開




時間発展のシュレディンガー方程式

微小時間後の波動関数

その展開係数としての飛び移り積分