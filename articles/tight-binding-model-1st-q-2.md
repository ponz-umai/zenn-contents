---
title: "第一量子化のTight-bindingモデル（後編）"
emoji: "🙌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: []
published: false
---


## Tight-bindingな近似その2：原子軌道が強く束縛されている感じの近似

さて、ここでより「強く束縛された状態」の描像に基づいた近似、つまり原子準位の広がりが小さいという近似を行っていきます。これにより簡単な場合であれば具体的なエネルギーバンド・固有関数が得られます。

### 重なり積分の近似



まず、先ほども考えたように、原子軌道$varphi_m(\boldsymbol{r}) $に対して、異なる格子点を中心に持つ原子軌道同士の「重なり」が小さいとして、思い切って次の積分を$0$と考えます。つまり


$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r} \simeq \delta_{m,n}\delta_{\boldsymbol{R},\boldsymbol{R}'}
$$

とします。この積分（Overlap Integral、重なり積分）が小さいというのは格子間距離に比べてWannier関数を近似している原子軌道の広がりが十分小さい、つまり「強く束縛されている」ようなイメージに対応しています。

これによりLCAO近似で表したWannier関数も正規直交性を持ち、Bloch関数も正規直交性を満たすようになります。

そして、先ほど得た永年方程式において、行列$S_{\boldsymbol{k}}$が単位行列になる：

$$
\begin{align*}
\left(S_{\boldsymbol{k}}\right)_{nm} &= \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
\\

&=

\sum_{\boldsymbol{R}'}
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}'-\boldsymbol{R})}\int
 \phi_n^*(\boldsymbol{r}-\boldsymbol{R}')
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}\\

 &=
 \sum_{\boldsymbol{R}'}
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}'-\boldsymbol{R})}
\delta_{n,m}\delta_{\boldsymbol{R}',\boldsymbol{R}}\\

&=
N\delta_{n,m}

\end{align*}
$$

ことで、永年方程式は

$$
\begin{Vmatrix}
(M_{\boldsymbol{k}})_{nm} - N\varepsilon_{n,\boldsymbol{k}}\delta_{nm}
\end{Vmatrix}
=0
$$

となり、行列$M_{\boldsymbol{k}}/N$の対角化に帰着します。

### 飛び移り積分

さらに

$$
\begin{align*}
\left(M_{\boldsymbol{k}}\right)_{nm} 
&=
 \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\hat{H}\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\

&=
\sum_{\boldsymbol{R}_i}
\sum_{\boldsymbol{R}_j}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}_i-\boldsymbol{R}_j)}\int
 \phi_n^*(\boldsymbol{r}-\boldsymbol{R}_i)
 \left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}_l}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\
\end{align*}
$$

とできますが、ここで並進の対称性を利用して$N$個の$\boldsymbol{R}_i$に対して積分は同じ値を取るので、$\boldsymbol{R}_i = \boldsymbol{0}$の場合の$N$倍を考えることで

$$
\begin{align*}
&=
N
\sum_{\boldsymbol{R}_j}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}_j)}\int
 \phi_n^*(\boldsymbol{r})
 \left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}_l}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\

\end{align*}
$$

となります。さらに原子軌道関数が孤立原子のハミルトニアンの固有関数であることを利用し、また$\boldsymbol{R}_l = \boldsymbol{R}_j$とそれ以外の場合で分けることで

$$
\begin{align*}

 &=
N
\sum_{\boldsymbol{R}_j}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}
\varepsilon_{n}^{\rm a}
\int
 \phi_n^*(\boldsymbol{r})
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}_j}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}

\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\

 &=
N
\varepsilon_{n}^{\rm a}
\delta_{n,m}\\

&\>\>\>\>+
N

\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r})d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}_j}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_j)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}_j}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}
\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}, \boldsymbol{R}_j}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}\\
\end{align*}
$$

となります。

ここで、第2項の積分を

$$
\begin{align*}
\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r})d\boldsymbol{r}
 
 &\equiv
\Delta\varepsilon_{nm}
 \end{align*}
$$

と定義します。この積分は（グロッソ・パラビチニによると）「結晶場積分」と呼ばれているそうなのですが、なんとググっても1件もヒットしないので、あまり一般的な用語ではないのかもしれません。（"Crystal Fiels Integral(s)"では1000件ほど（のみ）ヒット。）

なお、原点以外に中心を持つ局所ポテンシャルと、原点を中心とする原子軌道の内積の積を積分していますが、「中心からそれなりに離れた場所では原子ポテンシャルは大体一定になっているだろう」という近似のもと、

$$
\begin{align*}
\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}_l\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}_l)
   \right)
 \phi_m(\boldsymbol{r})d\boldsymbol{r}
 
 &\simeq
\Delta\varepsilon\delta_{nm}
 \end{align*}
$$

と定数の対角行列と近似されることもあるようです。このようにすればこの項はエネルギーを定数シフトさせるだけの効果となり、後の計算では無視していけるようになります。（本章では一応残しておきます。）

第3項は、原点と異なる格子点と、それと等しい局所ポテンシャルの積で、関数の局在性から近い位置の間の積分のみ値を持つと考え、特に原点から**再隣接格子点**へのベクトルを$\boldsymbol{R}_I$として

$$
\sum_{\boldsymbol{R}_j}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_j}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_j)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_j)d\boldsymbol{r}

 \simeq
\sum_{\boldsymbol{R}_I}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_I)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d\boldsymbol{r}

$$

と近似されます。ここで$\sum_{\boldsymbol{R}_I}$は、再隣接格子点のみ和を取ることを表します。

積分

$$
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_I)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d\boldsymbol{r}
 \equiv
 t_I^{n,m}
$$

は、上記のように$t$で表示されることが多く、 **"Transfer Integral"、"Hopping Integral"、「飛び移り積分」** 等と呼ばれます。（多分"transfer"の"t"だと思います）
多くの場合再隣接格子、時には次近接辺りまでの値が近似で採用されます。


また第4項は、全て中心が異なる関数の積分であり、同様にそれぞれの関数が局在していることから無視することができると考えます。（この形の積分を「3中心積分」等と呼ぶようです）

以上をまとめると、対角化すべき行列$M_{\boldsymbol{k}}/N$は

$$
\frac{M_{\boldsymbol{k}}}{N} =

\begin{bmatrix}
\varepsilon_{m_1}^{\rm a} - \Delta\varepsilon_{m_1m_1} & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_1m_2} -\Delta\varepsilon_{m_1m_2} & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_1m_3}- \Delta\varepsilon_{m_1m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_2m_1} - \Delta\varepsilon_{m_2m_1} &
  \varepsilon_{m_2}^{\rm a} -\Delta\varepsilon_{m_2m_2} &
   \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_2m_3} - \Delta\varepsilon_{m_2m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_3m_1} - \Delta\varepsilon_{m_3m_1} &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_3m_2} - \Delta\varepsilon_{m_3m_2} & 
 \varepsilon_{m_3}^{\rm a} - \Delta\varepsilon_{m_3m_3}
\end{bmatrix}
$$

となります。局所ポテンシャル$V(\boldsymbol{r})$や原子軌道関数$\phi_m(\boldsymbol{r})$は既に得られているものとしているので、行列の各要素は（少なくとも）数値的に得られて、後はサイズの小さい行列を対角化することで、展開に用いたエネルギー的に近い原子軌道の数（今回は3個）のエネルギーバンドと、3つの原子軌道の重ね合わせで表現される固有関数を得ることができます！

### 原子軌道や分子軌道の名前を用いたバンドの呼び方

この時、例えば展開に$s$軌道を選んだ場合、そこから得られるバンドを「$s$バンド」や「$s$的なバンド」、エネルギーが縮退した3つの$p$軌道で展開した場合に得られるバンドを「$p$バンド」「$p$的なバンド」などと呼びます。
また、（私はあまり詳しくないですが）展開の種類によっては分子軌道の用語「$\sigma$バンド」や「$\pi$バンド」なんかもよく見かけます。ここで「$\sigma$バンド」は$s, p_x, p_y$の展開、「$\pi$バンド」は残りの$p_z$軌道による展開に対応しているようです。


# Tight-binding Model

さて、上記のようにして固有値を求める方程式を得ることができました。この時、結晶場積分の値が全て等しいとして$\Delta\varepsilon\delta_{nm}$とおき、かつ展開に用いた軌道関数はエネルギーが近いものでしたので、思い切って原子準位$\varepsilon_{m_i}^{\rm a}$も全て等しいと近似し、その値を$\varepsilon^{\rm a}$と置けば、単位行列$I$を用いて

$$
\frac{M_{\boldsymbol{k}}}{N} = (\varepsilon^{\rm a} - \Delta\varepsilon)I + 
\begin{bmatrix}
0 & \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_1m_2}  &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_1m_3}\\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_2m_1} &
  0 &
   \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_2m_3} \\

\sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_3m_1}  &
 \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I^{m_3m_2}  & 
 0
\end{bmatrix}
$$

となり、エネルギーバンドの形状は飛び移り積分$t_I^{nm}$の値（と結晶構造）のみで決定されることになります。つまり飛び移り積分をパラメータとして見ることができるようになります。

そして、ここからは統計力学の範疇になるのですが、エネルギーバンドの形状が分かれば、様々な固体の物理量を計算することができるようになります。従ってエネルギーバンドは、固体の性質を端的に表現する関数であり、固体を表す「モデル」だともいえるわけで、上記のように考えるとそのモデルがパラメータとして$t$を用いて表現された形になってると見れます。

というわけで上記のような$t$をパラメータとして表したエネルギーバンドを指して、Tight-binding Modelということもあります。


## 簡単な例でのエネルギーバンド＝Tight-bindingモデル

最後に、最も簡単な例として、1次元格子において、固有関数を$s$軌道を用いて展開した場合について具体的に考えてみます。$s$軌道は縮退がないので、もっとも粗い近似では1つの原子軌道（Bloch和）のみで、Wannier関数（Bloch関数）を展開できます。

なお、1次元系をどのように表現すればよいかちょっとややこしいのですが、最も手っ取り早い方法は基本格子ベクトル$\boldsymbol{a}_i$を、

$$
|\boldsymbol{a}_2| = |\boldsymbol{a}_3| = A \gg |\boldsymbol{a}_1|
$$

と置いてしまうことかと思います。これは実際の1次元系について、それ以外の方向との格子間距離が離れている状況を表しています。（2方向に結晶サイズ1、つまり$N_2 = N_3 =1$の周期的境界条件を課しても上手くいきそうな気もするのですが、なんかややこしいのでやめておきます。）こうすると後述するように、考慮しない方向の飛び移り積分が$0$となり、1次元的なバンドが得られます。

さて、格子間隔が開いたとはいえ、固有関数（Bloch関数）は先述のように一つの$s$軌道関数を用いて、

$$
\varphi_{s,\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\phi_s(\boldsymbol{r}-\boldsymbol{R})
$$

と表せます。ここでややこしいですが、格子は1次元でも、その中の（固有）関数は3次元的に広がっており、展開に用いる原子軌道も3次元関数です。

1つの軌道で展開した場合は行列$M_{\boldsymbol{k}}$は行列ではなくなり、固有値方程式は


$$
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} + \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I
 - \varepsilon_{s,\boldsymbol{k}}
=0
\\

\Rightarrow
\varepsilon_{s,\boldsymbol{k}}
=
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} + \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_I
$$

と、すでに対角化された等式になります。
ここで$\varepsilon_{s}^{\rm a}$は$s$軌道の原子準位、$\Delta\varepsilon_{s}$は$s$軌道の結晶場積分、そして$\varepsilon_{s,\boldsymbol{k}}$が、電子の固有状態が持つ固有エネルギーの式（エネルギーバンド）を表します。


ここでここまではまだ3次元っぽいですが、ここで飛び移り積分$t_I$について考えると、1次元系っぽくなります。
すなわち、$\boldsymbol{a}_2, \boldsymbol{a}_3$方向の再隣接格子へのベクトルを$\boldsymbol{A}$とし、その時の飛び移り積分を$t_A$と置くと、

$$
t_A = \int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} \pm \boldsymbol{A})\phi_s(\boldsymbol{r}\pm \boldsymbol{A})d\boldsymbol{r}
$$

となりますが、ここで$A$は大きい数と置いたことから、$\boldsymbol{a}_2, \boldsymbol{a}_3$方向の飛び移り積分は$t_A\simeq0$となります。この結果エネルギーバンドは

$$
\varepsilon_{s,\boldsymbol{k}}
\simeq
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} + \sum_{R_1 = \pm a}e^{ik_1 R_1} t_I
$$

と、波数空間における1次元関数（$k_2, k_3$方向には定数関数）のエネルギーバンドが得られます。

さらに、$V(\boldsymbol{r})$が球対称ポテンシャルを仮定しており、また$s$軌道も角度部分は定数、つまり球対称であるため、$R_I = \pm a$いずれも同じ積分の値を持ちます。
最後に、$s$軌道は実数関数かつ正、つまり$\phi_s(\boldsymbol{r})^* = \phi_s(\boldsymbol{r}) > 0$であることと、$V$は引力ポテンシャルであることから$V(\boldsymbol{r})<0$より、積分は負の値になります。

そこで

$$
-\int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} + \boldsymbol{a}_1)\phi_s(\boldsymbol{r} + \boldsymbol{a}_1)
=

-\int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} - \boldsymbol{a}_1)\phi_s(\boldsymbol{r} - \boldsymbol{a}_1)
\equiv t >0
$$

と置いて、1次元関数の固有エネルギー（エネルギーバンド）は

$$
\begin{align*}
\varepsilon_{s}(k)
&=
\varepsilon_{s}^{\rm a} -\Delta\varepsilon_{s}-t \sum_{R_I = \pm a}e^{ikR} \\

&=
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} -2t \cos(ka)
\end{align*}
$$

と、$s$軌道の原子順位$\varepsilon_s^{\rm a}$（から結晶場積分を引いたもの）を中心に、飛び移り積分$t$程度に広がったバンド構造をとることがわかります。
