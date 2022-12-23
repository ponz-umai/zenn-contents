---
title: "固体・周期ポテンシャル内の電子状態～ブロッホの定理、エネルギーバンド、ブリルアンゾーン～"
emoji: "💠"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: false
---
# はじめに
前章で結晶の周期を持つ関数のFourier展開を定義することができたので、いよいよ固体、つまり周期的ポテンシャル内の電子の性質を考えていきます。

本稿のゴールであるTight-bindingモデルでは、格子点に孤立原子のポテンシャルが周期的に存在していると考えますが、本章ではより一般的に、ハミルトニアンが結晶の周期性を持つ場合にその固有状態が持つべき性質について整理します。

# 固体中の電子

## 固体内の電子のシュレーディンガー方程式

[以前の章](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)でも記載した通り、固体内の状態を色々と簡略化して考えると、原子核や内殻電子から受けるポテンシャルを全部まとめて1電子ポテンシャルに近似し、価電子間の相互作用を（一旦）無視することで、周期的に並んだ有効1電子ポテンシャルのもとを運動する相互作用しない多数の電子の状態を考えることになります。（下図も再掲しておきます）

![](/images/tb/solid.png)

多体のハミルトニアンは、位置$\boldsymbol{R}$の原子核＋内殻電子から、$i$番目の電子への有効ポテンシャルを$V_i(\boldsymbol{r}_i-\boldsymbol{R})$として、

$$
\mathcal{H}\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) = E\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots),\\
\mathcal{H} = \sum_i\left(  -\frac{\hbar^2}{2m}\nabla_i^2\right) + \sum_i\sum_{\boldsymbol{R}}V_i(\boldsymbol{r}_i-\boldsymbol{R}) 
$$

となります。

ここから変数分離をして、実際に扱うのは1電子シュレーディンガー方程式

$$
\left\{
    -\frac{\hbar^2}{2m}\nabla_i^2
    +
    \sum_{\boldsymbol{R}}V_i(\boldsymbol{r}_i-\boldsymbol{R})
\right\}
\varphi(\boldsymbol{r}) = \varepsilon\varphi(\boldsymbol{r})
$$

へと帰着されます。そこから求めた一体の固有関数をエネルギーが低い順からスレーター行列式に詰めて行けば基底状態の波動関数を得ることができます。

しかし1電子のシュレーディンガー方程式になったところで、やはり上式はそのままでは解けません。具体的な解を求めるためにはさらなる近似が必要となります。

## ブロッホの定理

近似にも色々ありますが、そのうちの一つに多電子原子の章で活用した変分法、すなわち波動関数の解になりそうな試行関数を仮定し、仮定した関数が解となるための条件を考えていく方法があります。

多電子原子について考えた際は、水素原子のシュレーディンガー方程式（中心力ポテンシャルの下でのシュレーディンガー方程式）が厳密に解けていたこともあり、その波動関数を出発点に置き、それなりに上手くいくことを見ました。それでは今回のような周期的なポテンシャルのもとでの試行関数はどのようなものを選べば良いでしょうか。

その指針となるのが、表題のブロッホの定理です。

ブロッホの定理は、格子ベクトルの並進に対して周期的なハミルトニアンの下で波動関数が従うべき条件についての定理です。

そのようなポテンシャルを$U(\boldsymbol{r})$と置きます。$U(\boldsymbol{r})$は前章・前々章で定義した格子ベクトル

$$
\begin{align*}
\boldsymbol{n} &= n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3\\
\end{align*}
$$

に対して

$$
U(\boldsymbol{r} + \boldsymbol{n}) = U(\boldsymbol{r})
$$

を満たします。また微分演算子$\nabla$は、座標の並進に対して不変です。^[なお、有効ポテンシャルの和$\sum_{\boldsymbol{R}}V_i(\boldsymbol{r}_i-\boldsymbol{R})$も格子ベクトルの並進に対して対称性を持ちます。$U(\boldsymbol{r}) = \sum_{\boldsymbol{R}}V_i(\boldsymbol{r}_i-\boldsymbol{R})$として並進後の$\boldsymbol{R} - \boldsymbol{n}=\boldsymbol{R}'$と置くと$U(\boldsymbol{r} + \boldsymbol{n})=\sum_{\boldsymbol{R}}V_i(\boldsymbol{r}_i + \boldsymbol{n}-\boldsymbol{R})=\sum_{\boldsymbol{R}'}{}^{'}V_i(\boldsymbol{r}_i -\boldsymbol{R}')$ここで$\boldsymbol{R}'$に関する和$\sum_{\boldsymbol{R}'}{}^{'}$の範囲は、並進したベクトル$\boldsymbol{n}$だけ$\boldsymbol{R}$の和からずれますが、先ほどの周期的境界条件から結局元の和と一致します。]

初めに結論から書くと、この時、シュレーディンガー方程式

$$

\hat{H}\varphi(\boldsymbol{r}) = \varepsilon\varphi(\boldsymbol{r}),\\

\hat{H} = 
    -\frac{\hbar^2}{2m}\nabla^2
    +
   U(\boldsymbol{r})
$$

を満たす波動関数は以下の条件を満たす、というのがブロッホの定理です。

::: message
固有値、固有関数はブロッホ波数$\boldsymbol{k}$, バンド指標$n$でラベルされ、
$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}),\\
\varepsilon_{n,\boldsymbol{k}}
$$
となる。^[ここで、波数$\boldsymbol{k}$と位置$\boldsymbol{r}$が同じ関数に現れていることが意味わからんと感じるかもしれません。私は長らくなんじゃこりゃと思いながら雰囲気で使っていました。これは「Fourier展開した際に波数$\boldsymbol{k}$に代表されるような波数成分のみ持つ位置$\boldsymbol{r}$の関数」というような意味です。Fourier係数$c_{\boldsymbol{k}}$とは別物です。紛らわしいですね。その辺は証明を追えばわかると思います。それに対応してブロッホ波数$\boldsymbol{k}$も、電子の運動量ではないことに注意してください。（とはいえ運動量っぽい性質に対応しており、「結晶運動量」等と呼ばれています。そのため$\boldsymbol{k}$の文字が使われているようです。ただ、本稿ではそのあたりまでは触れないと思います）]

この時の波動関数をブロッホ関数と呼び、ブロッホ関数$\varphi_{n,\boldsymbol{k}}$は格子ベクトル$\boldsymbol{n}$に対して

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r} + \boldsymbol{n}) = e^{i\boldsymbol{k}\cdot\boldsymbol{n}}
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) 
$$

を満たす。この条件は、結晶の周期性を持つ関数$u(\boldsymbol{r} + \boldsymbol{n}) = u(\boldsymbol{r})$を用いて

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})= e^{i\boldsymbol{k}\cdot\boldsymbol{r}}u(\boldsymbol{r})
$$

とも書ける。

ここで、固有値$\varepsilon_{n,\boldsymbol{k}}$、ブロッホ波数$\boldsymbol{k}$は、逆格子空間の特別な単位胞である（第1）ブリルアンゾーンに含まれる波数ベクトルとして定義され、波数$\boldsymbol{k}$に関して、逆格子ベクトル$\boldsymbol{K}$の周期関数となる。
:::

また、上記のブロッホの定理で導かれた固有値を考えることで、固体内の電子が持つエネルギーの性質を表すエネルギーバンドという概念が登場することになります。

本章ではそのあたりの内容について、例によって最低限、整理していきます。

以下、証明を書いていきます。

証明には並進の演算子を用いる方法と、波動関数が満たす条件を具体的に求める方法の2パターンあり、本章では後者の方針で証明を記載します。（前者の方がすっきりと書けるのですがその前提となる定理を色々と追記する必要があり、地道に理解できる後者の方向で書いていきます）

## ブロッホの定理の前半部分の証明（1次元）

初めに、色々と考えやすくなるので1次元の場合について証明を記載します。3次元の場合も後に触れますが、ほぼ同じです。

まずは「ブロッホ波数」とかいうものについては置いておいて、周期的なハミルトニアンの固有関数が格子ベクトル$\boldsymbol{n}$に対して

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r} + \boldsymbol{n}) = e^{i\boldsymbol{k}\cdot\boldsymbol{n}}
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) 
$$

または、結晶の周期性を持つ関数$u(\boldsymbol{r} + \boldsymbol{n}) = u(\boldsymbol{r})$を用いて

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})= e^{i\boldsymbol{k}\cdot\boldsymbol{r}}u(\boldsymbol{r})
$$

を満たすことを確かめていきます。

### 波数表示でのシュレーディンガー方程式

1次元周期ポテンシャルの下でのシュレーディンガー方程式を以下のように書きましょう。

$$
\left(
    \frac{-\hbar^2}{2m}\frac{d^2}{d x^2} + U(x)
\right)
\varphi(x) = \varepsilon\varphi(x)
$$

格子間隔（格子定数）を$a$、$n$を整数として、ポテンシャルは周期性

$$
U(x + na) = U(x)
$$

を満たします。
ここで境界条件として、十分大きな数$N$を用いて

$$
\varphi(x + Na) = \varphi(x)
$$

を考えます。^[この条件はボルンーフォンーカルマンの周期的境界条件と呼ばれているようです。]

この時、ポテンシャル$U(x)$は

$$
U(x) = \sum_K U_K e^{iKx},\\
K = \frac{2\pi}{a}l, \>\>l = 0,\pm 1, \pm 2\cdots
$$

とFourier展開できます。ここで波数$K$は逆格子に対応するものなので、その意味を強調するために大文字で書いています。

また、周期的境界条件を設定したことから波動関数も周期$Na$を持ち、

$$
\varphi(x) = \sum_q c_q e^{iqx},\\
q = \frac{2\pi}{Na}m, \>\>m = 0,\pm 1, \pm 2\cdots
$$

と展開できます。

これらをシュレーディンガー方程式に代入して整理します。運動エネルギー項に対しては

$$
\frac{d^2}{d x^2}\sum_{q}c_{q}e^{iqx}
=
\sum_{q}q^2c_{q}e^{iqx}
$$

となります。またポテンシャル項について、

$$
\begin{align*}
U(x)\varphi(x) 

&= 
\sum_{K}\sum_{q}U_{K}
c_{q}
e^{i(K+q)x}\\

&=
\sum_{K}\sum_{q'}U_{K}
c_{q'-K}
e^{iq'x}\\
\end{align*}
$$

ですが、$q + K = \frac{2\pi}{Na}(m + lN)\rightarrow q'$の置き換えに対しても$\sum_{q'}$の和の範囲は変わらないので、改めて$q'\rightarrow q$と置きなおして、

$$
U(x)\varphi(x) 
=
\sum_{K}\sum_{q}U_{K}
c_{q-K}
e^{iqx}\\
$$

と書けます。これらをまとめて、波数表示のシュレーディンガー方程式を

$$
\sum_{q}e^{iqx}
\left[
    (\frac{\hbar^2}{2m}q^2 - \varepsilon)c_{q} + 
    \sum_{K}U_{K}
c_{q-K}

    \right]
    =0
$$

と書き換えることができ、Fourier係数$c_{q}$の満たす方程式

$$
 (\frac{\hbar^2}{2m}q^2 - \varepsilon)c_{q} + 
    \sum_{K}U_{K}
c_{q-K} = 0
$$

を得ることができます。

### ブロッホ関数の導出

この方程式は任意に選んだある$q$と、その逆格子分異なる$q-K, q\pm K', q\pm K''$の係数を持つFourier係数のみを含む連立方程式のブロックになっています。（このようなブロックを代表する$q$は、逆格子の周期に関して最小単位となる逆格子空間の単位胞内でとればよいことを後に確かめます）

行列表示の都合上、逆格子$K=\frac{2\pi}{a}m$は、$m=0,\pm 1, \pm 2\cdots$なので、それぞれ$K_0, \pm K_1, \cdots,\pm K_m\cdots$と書くと、ある$q$を選ぶと、その逆格子分だけ異なる係数に関する行列の固有値方程式

$$
\begin{bmatrix}
\ddots &   &  &  &   \\
 & \ddots & \\
  & & \ddots  & & &\\
 & & & \ddots \\
\cdots & \cdots & U_{-K_2} &  U_{-K_1} & \frac{\hbar^2}{2m}q^2 + U_0 & U_{K_1} & U_{K_2}&\cdots \\
\\
   & & &&& \ddots & \\
&& & & &&\ddots &  \\
&&& & & &&\ddots &  \\
\end{bmatrix}

\begin{bmatrix}
\vdots \\
c_{q - K_m} \\
\vdots \\
c_q \\
\vdots \\
c_{q + K_m} \\
\vdots
\end{bmatrix}

= \varepsilon 
\begin{bmatrix}
\vdots \\
c_{q - K_m} \\
\vdots \\
c_q \\
\vdots \\
c_{q + K_m} \\
\vdots
\end{bmatrix}

$$

がブロックで抜き出されることになります。この方程式は行列$H_q$を

$$
H_q = \begin{bmatrix}
\ddots &   &  &  &   \\
 & \ddots & \\
  & & \ddots  & & &\\
 & & & \ddots \\
\cdots & \cdots & U_{-K_2} &  U_{-K_1} & \frac{\hbar^2}{2m}q^2 + U_0 & U_{K_1} & U_{K_2}&\cdots \\
\\
   & & &&& \ddots & \\
&& & & &&\ddots &  \\
&&& & & &&\ddots &  \\
\end{bmatrix}
$$

として固有値方程式

$$
|H_q -\varepsilon I| = 0
$$

より無限個の固有値$\varepsilon_{n,q}$を求めれば、その固有値に対応した固有ベクトル

$$
H_q
\begin{bmatrix}
\vdots \\
c_{n, q - K_m} \\
\vdots \\
c_{n,q} \\
\vdots \\
c_{n, q + K_m} \\
\vdots
\end{bmatrix}
=
\varepsilon_{n,q}
\begin{bmatrix}
\vdots \\
c_{n, q - K_m} \\
\vdots \\
c_{n,q} \\
\vdots \\
c_{n, q + K_m} \\
\vdots
\end{bmatrix}
$$

を得ることができます。


この結果から、元々求めたかったシュレーディンガー方程式に対する固有関数を先ほど求めた固有ベクトルの係数$\{c_{n,q-Km}\}$を用いて

$$
\varphi_{n,q}(x) = \sum_{m}c_{n,q-K_m}e^{i(q-K_m)x}
$$

を得ます。
すなわち、周期的なハミルトニアンに対する固有値は、$n$でラベルされた様々な波数$q-K_m$の係数の組み$c_{n,q-K_m}$のみをFourier成分として持つ関数として表される、ということがわかります。


実際、$\varphi_{n,q}$は$\hat{H} = \frac{-\hbar^2}{2m}\frac{d^2}{d x^2} + U(x)$に対して

$$
\begin{align*}
\hat{H}\sum_m c_{n,q-K_m}e^{i(q-K_m)x} 
&= 
\sum_m\frac{\hbar^2}{2m}(q-K_m)^2c_{n,q-K_m}e^{i(q-K_m)x} + \sum_m \sum_{m'} U_{K_{m'}}c_{n,q-K_m}e^{i(q-(K_m - K_{m'}))x}\\

&=
\sum_m\frac{\hbar^2}{2m}(q-K_m)^2c_{n,q-K_m}e^{i(q-K_m)x}  + \sum_{K_{m''}} \sum_{K_{m'}} U_{K_{m'}}c_{n,q-K_{m''}-K_{m'}}e^{i(q- K_{m''})x}\\

&=
\sum_m\frac{\hbar^2}{2m}(q-K_m)^2c_{n,q-K_m}e^{i(q-K_m)x}  + \sum_{K_m} \sum_{K_{m'}} U_{K_{m'}}c_{n,q-K_m-K_{m'}}e^{i(q- K_m)x}\\

&=
\sum_m e^{i(q-K_m)x} \left[
    \frac{\hbar^2}{2m}(q-K_m)^2c_{n,q-K_m} + \sum_{m'} U_{K_{m'}}c_{n,q-K_m-K_{m'}}
    \right]
    \\

&=
\sum_m e^{i(q-K_m)x} \left[
    \sum_l(H_q)_{m,l} 
    c_{n,q-K_l} 
    \right]\\

&=
\sum_m e^{i(q-K_m)x} \left[
    \varepsilon_{n,q} 
    c_{n,q-K_m} 
    \right]\\

&=
\varepsilon_{n,q}\varphi_{n,q}(\boldsymbol{r})

\end{align*}
$$

と、固有状態になっていることがわかります。

またラベル$m$についての和は、逆格子$K$についての和と等しいので、エネルギー固有値$\varepsilon_{n,q}$に対応する固有関数は

$$
\varphi_{n,q}(x) = \sum_K c_{n,q-K}e^{i(q-K)x}
$$

とも書けます。



### ブロッホ関数の座標並進に対する性質（ブロッホの定理の前半の証明）

さて、当初目標のブロッホの定理の前半部分ですが、これは実際に代入してみれば確かめることができます。

この関数は格子の並進$x\rightarrow x + na$に対して、$e^{iKna} = 1$より、

$$
\begin{align*}
    \varphi_{n,q}(x + na) &= \sum_K c_{n,q-K}e^{i(q-K)x}e^{i(q-K)na}\\
    &=
    \sum_K c_{n,q-K}e^{i(q-K)x}e^{iqna}\\
    &=
    e^{iqna}\varphi_{n,q}(x) 
\end{align*}
$$

を満たします。また、

$$
\varphi_{n,q}(x) = e^{iqx}\sum_K c_{n,q-K}e^{-iKx}
$$
と書くと、右辺$e^{iqx}$を除いた部分は

$$
\sum_K C_{n,q-K}e^{-iK(x + na)} = \sum_K C_{n,q-K}e^{-iKx}e^{-iKna}
=
\sum_K C_{n,q-K}e^{-iK(x)}
$$

を満たし、結晶の周期性を持つ関数$u(x)$と書けます。

以上で、ブロッホの定理として示した前半部分

::: message
ブロッホ関数$\varphi_{n,\boldsymbol{k}}$は格子ベクトル$\boldsymbol{n}$に対して

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r} + \boldsymbol{n}) = e^{i\boldsymbol{k}\cdot\boldsymbol{n}}
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) 
$$

を満たす。この条件は、結晶の周期性を持つ関数$u(\boldsymbol{r} + \boldsymbol{n}) = u(\boldsymbol{r})$を用いて

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})= e^{i\boldsymbol{k}\cdot\boldsymbol{r}}u(\boldsymbol{r})
$$
とも書ける
:::

を示すことができました。「ブロッホ波数」についての説明や波数空間での周期性の話に移る前に、先に3次元版でも同様に考えておきます。


## 3次元の場合のブロッホの定理の証明 

上記で1次元の場合についてブロッホの定理を示しましたが、3次元の場合も以下の点を変えれば同じ流れで考えることができます。

ポテンシャル$U(\boldsymbol{r})$のFourier展開は前章で行ったように逆格子ベクトル$\boldsymbol{K}$：

$$
    \boldsymbol{K} 

=
 m_1\frac{2\pi\boldsymbol{a}_2\times\boldsymbol{a}_3}{\boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)}
+

 m_2\frac{2\pi\boldsymbol{a}_3\times\boldsymbol{a}_1}{\boldsymbol{a}_2\cdot(\boldsymbol{a}_3\times\boldsymbol{a}_1)}
+
 m_3\frac{2\pi\boldsymbol{a}_1\times\boldsymbol{a}_2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}
$$

を用いた

$$
U(\boldsymbol{r}) = \sum_{\boldsymbol{K}}U_{\boldsymbol{K}}e^{i\boldsymbol{K}\cdot\boldsymbol{r}}
$$

となります。

また、波動関数$\varphi(\boldsymbol{r})$の境界条件として、先ほど考えた条件の3次元版、すなわち基本併進ベクトルを$\boldsymbol{a}_1,\boldsymbol{a}_2,\boldsymbol{a}_3$としたとき、十分大きな数$N_1, N_2, N_3$ （$N_1\times N_2\times N_3 \times v_c$ (単位胞の体積）が大体結晶のサイズになるくらい）を考え、

$$
\varphi(\boldsymbol{r} + N_i\boldsymbol{a}_i) = \varphi(\boldsymbol{r})
$$

という条件を考えます。すると「逆格子ベクトル」^[かっこをつけている理由は、正確には「逆格子ベクトル」は結晶の格子ベクトルの双対ベクトルに対して使用する用語だと教科書に書いてあったからです。なので境界条件の周期性に対して「逆格子ベクトル」と書くのは正確ではないっぽいのですが、ここではわかりやすさを優先してこのように書きます]

$$
\begin{align*}
    \boldsymbol{q} 

&=
 m_1\frac{2\pi N_2\boldsymbol{a}_2\times N_3\boldsymbol{a}_3}{ N_1\boldsymbol{a}_1\cdot( N_2\boldsymbol{a}_2\times N_3\boldsymbol{a}_3)}
+

 m_2\frac{2\pi N_3\boldsymbol{a}_3\times N_1\boldsymbol{a}_1}{ N_2\boldsymbol{a}_2\cdot( N_3\boldsymbol{a}_3\times N_1\boldsymbol{a}_1)}
+
 m_3\frac{2\pi N_1\boldsymbol{a}_1\times N_2\boldsymbol{a}_2}{ N_3\boldsymbol{a}_3\cdot( N_1\boldsymbol{a}_1\times N_2\boldsymbol{a}_2)}\\

 &=
 
 \frac{m_1}{N_1}\frac{2\pi\boldsymbol{a}_2\times\boldsymbol{a}_3}{\boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)}
+

\frac{m_2}{N_2}\frac{2\pi\boldsymbol{a}_3\times\boldsymbol{a}_1}{\boldsymbol{a}_2\cdot(\boldsymbol{a}_3\times\boldsymbol{a}_1)}
+
\frac{m_3}{N_3}\frac{2\pi\boldsymbol{a}_1\times\boldsymbol{a}_2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}

\end{align*}
$$

の平面波$e^{i\boldsymbol{q}\cdot \boldsymbol{r}}$を用いて以下のように展開できます：

$$
\varphi(\boldsymbol{r}) = \sum_{\boldsymbol{q}}c_{\boldsymbol{q}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}
$$

これらをシュレーディンガー方程式に代入して整理すると、運動エネルギー項に対しては

$$
\nabla^2\sum_{\boldsymbol{q}}c_{\boldsymbol{q}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}
=
\sum_{\boldsymbol{q}}q^2c_{\boldsymbol{q}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}
$$

となります。またポテンシャル項について、

$$
\begin{align*}
U(\boldsymbol{r})\varphi(\boldsymbol{r}) 

&= 
\sum_{\boldsymbol{K}}\sum_{\boldsymbol{q}}U_{\boldsymbol{K}}
c_{\boldsymbol{q}}
e^{i(\boldsymbol{K}+\boldsymbol{q})\cdot\boldsymbol{r}}\\

&=
\sum_{\boldsymbol{K}}\sum_{\boldsymbol{q}'}U_{\boldsymbol{K}}
c_{\boldsymbol{q}'-\boldsymbol{K}}
e^{i\boldsymbol{q}'\cdot\boldsymbol{r}}\\
\end{align*}
$$

ですが、$\boldsymbol{q}+\boldsymbol{K}\rightarrow\boldsymbol{q}'$の置き換えに対しても$\sum_{\boldsymbol{q}'}$の和の範囲は変わらないので、改めて$\boldsymbol{q}'\rightarrow\boldsymbol{q}$と置きなおして、

$$
U(\boldsymbol{r})\varphi(\boldsymbol{r}) 
=
\sum_{\boldsymbol{K}}\sum_{\boldsymbol{q}}U_{\boldsymbol{K}}
c_{\boldsymbol{q}-\boldsymbol{K}}
e^{i\boldsymbol{q}\cdot\boldsymbol{r}}\\
$$

と書けます。これらをまとめて、波数表示のシュレーディンガー方程式を

$$
\sum_{\boldsymbol{q}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}
\left[
    (\frac{\hbar^2}{2m}q^2 - \varepsilon)c_{\boldsymbol{q}} + 
    \sum_{\boldsymbol{K}}U_{\boldsymbol{K}}
c_{\boldsymbol{q}-\boldsymbol{K}}

    \right]
    =0
$$

と書き換えることができ、Fourier係数$c_{\boldsymbol{q}}$の満たす方程式

$$
 (\frac{\hbar^2}{2m}q^2 - \varepsilon)c_{\boldsymbol{q}} + 
    \sum_{\boldsymbol{K}}U_{\boldsymbol{K}}
c_{\boldsymbol{q}-\boldsymbol{K}} = 0
$$

を得ることができます。ここから、逆格子ベクトルに何らかの順番を考えて行列表示で考えれば1次元の場合と同様に、ある波数$\boldsymbol{q}$と、その逆格子ベクトルだけ異なるFourier級数からなるブロッホ関数

$$
\varphi_{n,\boldsymbol{q}}(\boldsymbol{r}) = \sum_{\boldsymbol{K}}c_{n,\boldsymbol{q} - \boldsymbol{K}}
e^{i(\boldsymbol{q} - \boldsymbol{K} )\cdot \boldsymbol{r}}
$$

を固有関数として導くことができます。

この関数に対して格子ベクトルの並進を考えることで、1次元の場合と同様に


$$
\varphi_{n,\boldsymbol{q}}(\boldsymbol{r} + \boldsymbol{n}) = e^{i\boldsymbol{k}\cdot\boldsymbol{n}}
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) 
$$

を満たすことが確認できます。さらに、

$$
u(\boldsymbol{r}) \equiv
\sum_{\boldsymbol{K}}c_{n,\boldsymbol{q} - \boldsymbol{K}}
e^{-i \boldsymbol{K} \cdot \boldsymbol{r}}
$$

が結晶の周期性を持つ

$$
u(\boldsymbol{r} + \boldsymbol{n}) = u(\boldsymbol{r})
$$

ことが示せます。



# ブリルアンゾーン

さて、ここまでラベル$q$または$\boldsymbol{q}$の選び方に関して特に注意を払っておりませんでしたが、以下のように考えることで、固有エネルギー・固有関数が波数$\boldsymbol{q}$に関して、逆格子ベクトルの周期性を持っていると考えることができます。

すなわち、先ほどブロッホの定理として示した内容の後半部分、

::: message
固有値、固有関数はブロッホ波数$\boldsymbol{k}$, 離散的なバンド指標$n$でラベルされ、
$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}),\\
\varepsilon_{n,\boldsymbol{k}}
$$
となる（固有値$\varepsilon_{n,\boldsymbol{k}}$は$\boldsymbol{k}$についての連続関数であることを表して$\varepsilon_n(\boldsymbol{k})$とも書かれる）。
ブロッホ波数$\boldsymbol{k}$は、逆格子空間の特別な単位胞である（第1）ブリルアンゾーンに含まれる波数ベクトルとして定義される。
:::

が成立していることを本節で確認し、その流れで「ブリルアンゾーン」という概念を導入します。

### 固有値・固有関数の波数空間での周期性（1次元版）（ブロッホの定理後半の証明）

初めに、例によってあまりややこしくない1次元の場合から考えていきます。先ほど考えた1次元の場合のブロッホ関数の導出を思い出してみると、あるラベル$q$を選ぶとその逆格子分だけ異なる波数のみから係数$c_{q+K}$に関する方程式ブロック

$$
\begin{bmatrix}
\ddots &   &  &  &   \\
 & \ddots & \\
  & & \ddots  & & &\\
 & & & \ddots \\
\cdots & \cdots & U_{-K_2} &  U_{-K_1} & \frac{\hbar^2}{2m}q^2 + U_0 & U_{K_1} & U_{K_2}&\cdots \\
\\
   & & &&& \ddots & \\
&& & & &&\ddots &  \\
&&& & & &&\ddots &  \\
\end{bmatrix}

\begin{bmatrix}
\vdots \\
c_{q - K_m} \\
\vdots \\
c_q \\
\vdots \\
c_{q + K_m} \\
\vdots
\end{bmatrix}

= \varepsilon 
\begin{bmatrix}
\vdots \\
c_{q - K_m} \\
\vdots \\
c_q \\
\vdots \\
c_{q + K_m} \\
\vdots
\end{bmatrix}

$$


が形成されるのでした。

一方、例えば先ほど選んだ$q$に対して、ある逆格子分だけ異なる$q' = q - K_l$に対する方程式

$$
 (\frac{\hbar^2}{2m}(q-K_l)^2 - \varepsilon)c_{q-K_l} + 
    \sum_{K'}U_{K'}
c_{q-K_l-K'} = 0
$$


を考えてみたとします。
そうした場合も、先ほどと同じように$q'$から全ての逆格子分だけ異なる係数からなる連立方程式ブロック

$$
\begin{bmatrix}
\ddots &   &  &  &   \\
 & \ddots & \\
  & & \ddots  & & &\\
 & & & \ddots \\
\cdots & \cdots & U_{-K_{l+2}} &  U_{-K_{l+1}} & \frac{\hbar^2}{2m}(q-K_l)^2 + U_{K_{l}} & U_{K_{l+1}} & U_{K_{l+2}}&\cdots \\
\\
   & & &&& \ddots & \\
&& & & &&\ddots &  \\
&&& & & &&\ddots &  \\
\end{bmatrix}

\begin{bmatrix}
\vdots \\
c_{q - K_l- K_m} \\
\vdots \\
c_{q-K_l} \\
\vdots \\
c_{q  - K_l + K_m} \\
\vdots
\end{bmatrix}

= \varepsilon 
\begin{bmatrix}
\vdots \\
c_{q - K_l- K_m} \\
\vdots \\
c_{q-K_l} \\
\vdots \\
c_{q  - K_l + K_m} \\
\vdots
\end{bmatrix}
$$

を得ることができますが、その中には必ず先に選んだ$q = q' + K_l$も含まれており、結局全ての方程式を参照すると全く同じ式の集まりであることがわかります。

当然、その方程式から作られる行列の固有値$\varepsilon_{n,q'}$も、固有関数$\varphi_{n,q'}(x)$も先に考えた固有値$\varepsilon_{n,q}$も、固有関数$\varphi_{n,q}(x)$と全く同じものになります。

これはある波数$q$を選んで得た（無限個の$n$でラベルされた）固有値$\varepsilon_{n,q}$、固有関数$\varphi_{n,q}$と、先ほどの波数から逆格子だけ異なる波数$q-K$を選んで得た固有値$\varepsilon_{n,q-K}$、固有関数$\varphi_{n,q-K}$が一致していることを意味しています。ラベル$q$は方程式を指定するラベルであり、逆格子分ずれたラベルからは同じ方程式が指定されるわけです^[ショッピングモールの入り口がたくさんあったとしても、どこから入っても店舗さえ指定されていれば同じ商品が売ってる、みたいなもんです。]。

そこで、ラベル$q$に関しては「逆格子の周期の最小単位の領域に含まれる$q$のみを考えればよい」ということになります。

図式的に表すと、下図のようにある波数、特に$k = 0$を選ぶと、それに対応して逆格子$\pm K_1, \pm K_2\cdots$だけ異なる波数$k = \pm K_1, \pm K_2\cdots$の方程式の組が得られるので、$k = 0$を中心に逆格子の周期性の最小単位の$k$を選べば全ての$k$に対する方程式が得られる、と考えられます。

![](/images/tb/1d-bz.png)

あるいは少し考えを拡張して、波数ベクトル$q$をラベルとする固有値、固有関数は、波数空間内で逆格子の周期性

$$
\varepsilon_{n,q + K} = \varepsilon_{n,q},\\
\varphi_{n,q + K}(x) = \varphi_{n,q}(x)
$$

を持つ周期関数であるとも考えることができます。

周期関数であれば、関数を考える際やグラフで表示する際、1周期分のみを考えれば良いわけです。

いずれにしても、波数空間の周期の最小単位内の波数を代表して選べば、全ての波数に対する方程式が得られることがわかります。
そのような波数空間における単位胞内の波数を「ブロッホ波数」と呼び、記号$k$で表します。


そして周期の最小単位の領域のことをブリルアンゾーンと呼び、特に原点を中心として、隣接する逆格子点との垂直二等分面（2次元だと線、1次元だと点）で区切られた領域を第一ブリルアンゾーンと呼びます。


### 固有値・固有関数の波数空間での周期性（3次元版）

3次元で考えた場合も、先ほどと同様に、ある波数ベクトル$\boldsymbol{q}$に対して得られる方程式ブロックと、逆格子ベクトルだけ異なる$\boldsymbol{q} - \boldsymbol{K}$に対して得られる方程式ブロックが等しく、その結果固有エネルギー$\varepsilon_{n,\boldsymbol{q}}$、ブロッホ関数$\varphi_{n,\boldsymbol{q}}(\boldsymbol{r})$が逆格子ベクトルの周期を持つことがわかります。

さて、1次元の場合はただ隣接する逆格子点までの中点にはさまれた領域を考えればよかったですが、3次元の場合はどのような領域を考えればよいでしょうか。
素朴には、最小単位であれば良いので逆格子ベクトル

$$
\begin{align*}
    \boldsymbol{K} 

&=
 m_1\frac{2\pi\boldsymbol{a}_2\times\boldsymbol{a}_3}{\boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)}
+

 m_2\frac{2\pi\boldsymbol{a}_3\times\boldsymbol{a}_1}{\boldsymbol{a}_2\cdot(\boldsymbol{a}_3\times\boldsymbol{a}_1)}
+
 m_3\frac{2\pi\boldsymbol{a}_1\times\boldsymbol{a}_2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}\\

&=m_1\boldsymbol{b}_1 + m_2\boldsymbol{b}_2 + m_3\boldsymbol{b}_3
\end{align*}
$$

に対して、領域

$$
\boldsymbol{k} = s\boldsymbol{b}_1 + t\boldsymbol{b}_2 + u\boldsymbol{b}_3,\\
0\leq s < 1, 0\leq t < 1, 0\leq u < 1
$$

の単位胞内の$\boldsymbol{k}$を考えればあとは逆格子ベクトルを足し引きすれば全ての波数$\boldsymbol{q}$をカバーできるわけですが、それよりももっと性質の良い単位胞の取り方として、1次元版でちらっと述べたように
「原点から隣接する全ての逆格子ベクトルに対する垂直2等分面を考え、その面に囲まれる領域」
という単位胞を考えます。このような単位胞を「Wigner-Seitz cell（ウィグナー・サイツ・セル）」と呼び、特に波数空間におけるWigner-Seitz cellをブリルアンゾーン、その原点から最も近い領域を第一ブリルアンゾーンと呼びます。

そして、このような（第一）ブリルアンゾーン内で定義された波数ベクトルを「ブロッホ波数」と呼び、記号$\boldsymbol{k}$を用いて表します。



ここはあまり深く突っ込まないのですが、なぜこんなややこしいことをするかというと、Wigner-Seitz cellはあるブラベー格子（これは実空間の格子であっても逆格子空間の逆格子ベクトルの集まりであっても）を考えた際に、そのブラベー格子が持つ対称性と等しい対称性を持つ、という性質があるからのようです。

この辺はいずれ、ちゃんと勉強して何かしら追記しますが、この段階ではひとまずこんなところで。

ここまでの流れでわかるように（ブロッホ）波数は、「ある「連立方程式ブロック」を代表する波数」という意味ですので、ラベルが直接関数を指定していた自由電子の波動関数$\varphi_k(x) = Ae^{ikx}$や中心力場内の電子の固有関数$\varphi_{nml} = R_{nl}Y_l^m$とは若干違う意味合いで使われおり、「方程式群を指定するラベル」$k$と、「そこから得られる関数を指定するラベル$n$」というセットとなっております。
その上に自由電子と全く同じラベルが使われているので中々ややこしいなと感じますが、「結晶内の波数っぽい物理量に対応している」的な意味で自由電子の波数と同じ記号が使われ、「結晶運動量」等と呼ばれたりしているようです。（先述の通り本稿では今のところはここも深く突っ込まないことにします）

さて、話を戻して、以上で格子の周期を持つハミルトニアンの下での固有値、固有ベクトルの性質であるブロッホの定理

::: message
固有値、固有関数はブロッホ波数$\boldsymbol{k}$, バンド指標$n$でラベルされ、
$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}),\\
\varepsilon_{n,\boldsymbol{k}}
$$
となる。

この時の波動関数をブロッホ関数と呼び、ブロッホ関数$\varphi_{n,\boldsymbol{k}}$は格子ベクトル$\boldsymbol{n}$に対して

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r} + \boldsymbol{n}) = e^{i\boldsymbol{k}\cdot\boldsymbol{n}}
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) 
$$

を満たす。この条件は、結晶の周期性を持つ関数$u(\boldsymbol{r} + \boldsymbol{n}) = u(\boldsymbol{r})$を用いて

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})= e^{i\boldsymbol{k}\cdot\boldsymbol{r}}u(\boldsymbol{r})
$$

とも書ける。

ここで、固有値$\varepsilon_{n,\boldsymbol{k}}$、ブロッホ波数$\boldsymbol{k}$は、逆格子空間の特別な単位胞である（第1）ブリルアンゾーンに含まれる波数ベクトルとして定義され、波数$\boldsymbol{k}$に関して、逆格子ベクトル$\boldsymbol{K}$の周期関数となる。
:::

が成り立つことを確認できました。

# エネルギーバンド

最後に、エネルギー固有値$\varepsilon_{n,\boldsymbol{k}}$の持つ性質から、ブロッホ波数に対してプロットしたブロッホ電子の固有エネルギーの性質を確認して終わりにします。

表題の通り、エネルギー固有値$\varepsilon_{n,\boldsymbol{k}}$は、ブロッホ波数に対して「帯（バンド）」のようなグラフでプロットされます。

これはどのようにして$\varepsilon_{n,\boldsymbol{k}}$を求めたかを考えればよく、ブロッホ波数$\boldsymbol{k}$（第一ブリルアンゾーン内の波数に限ることを反映して、ここでは$\boldsymbol{k}$を使います）をパラメータとして持つ行列$H_{\boldsymbol{k}}$の固有値方程式

$$
|H_{\boldsymbol{k}} - \varepsilon I| = 0
$$

の（無限個の）根として$\varepsilon_{n,\boldsymbol{k}}$が得られたのでした。

この時、パラメータに依存した行列の固有値はそのパラメータに対して連続かつ微分可能である、らしいです。
**すいません。**このあたりの証明がよくわかっていないので、ここは曖昧なままにしてしまいます。いずれ線形代数か量子力学の最初の方か色々読み返して、ちゃんとやります。

またこちらも全く曖昧なのですが、行列$H_{\boldsymbol{k}}$の対角要素$\frac{\hbar^2}{2m}(\boldsymbol{q} - \boldsymbol{K})^2$の逆格子ベクトル$\boldsymbol{K}$が離散的であることから、対角化して得られる固有値$\varepsilon_n$は$n$に対して離散的になる、らしいです。

全く正確ではないと思うのですが、今のところは、行列式

$$
|H_{\boldsymbol{k}} - \varepsilon I|=0
$$

を$\varepsilon$の多項式として考えてみると、

$$
\cdots + A(\boldsymbol{k})\varepsilon^n + A'(\boldsymbol{k})\varepsilon^{n-1} + \cdots = 0
$$

みたいな方程式になるわけですが、この根は上記方程式左辺を$\varepsilon$の関数としてプロットした場合に横軸との交点になるわけで、まあ確かに、$\boldsymbol{k}$の変化に対して滑らかに交点は変化しそうですし、交点は離散的になりそうな気がします。

# おわりに


