---
title: "第一量子化のTight-bindingモデル（後編）"
emoji: "🧶"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# はじめに

前章で、固体内の周期（多体）ポテンシャルを周期的に並んだ1体ポテンシャル（孤立原子ポテンシャル）に近似し、その固有状態を用いて固体の固有関数を展開しました。

また、「強く束縛された状態」な仮定のもと、固有関数のうちあるバンドに属するBloch関数はお互いに近い原子準位
つまり$\varepsilon_{n_q,\boldsymbol{k}} \simeq \varepsilon_{m_1}^{\rm a}, \varepsilon_{m_2}^{\rm a}, \varepsilon_{m_3}^{\rm a}\cdots \varepsilon_{m_p}^{\rm a}$を持つ少数の原子軌道関数のみを用いて近似的に展開

$$
\varphi_{n_q,\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{m = m_1,m_2,\cdots m_p}b_m^{n_q}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

でき、その展開係数と固有エネルギーは以下のようなサイズの小さい一般化固有値方程式を解くことで得られることを見ました。


$$
\begin{vmatrix}
  M_{\boldsymbol{k}}^q - \varepsilon_{\boldsymbol{k}} S_{\boldsymbol{k}}^q
\end{vmatrix} = 0
$$


$$
\left(M_{\boldsymbol{k}}^q\right)_{ij} \equiv \int \Phi_{m_i,\boldsymbol{k}}^{q*}(\boldsymbol{r})\hat{H}\Phi^{q}_{m_j,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\
\left(S_{\boldsymbol{k}}^q\right)_{ij} \equiv \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

$$
\Phi_{m,\boldsymbol{k}}^q\equiv \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\phi_m(\boldsymbol{r} - \boldsymbol{R}), \>\>\>\>m = m_1, m_2, \cdots m_p
$$


ここまででも数値的に方程式を解けばエネルギーバンドは求められますし、あるいは3軌道くらいでの展開であれば3次方程式を解けばいいので力業で回析的な答えもわかるかもしれませんが、いかんせんそれでは何がどうなっているのかイメージがつきづらいですよね。

また、Tight-bindingモデルは単体で用いられるというよりも、Habbard モデル等さらに進んだモデルの出発点となるものです。しかしその解が数値的に得られないのでは、いささか使いにくさもあるのではないでしょうか。
従って解ければいいわけではなく、使いやすい形でモデル化したいわけです。

というわけで、より物理的なイメージをはっきりさせるために、あるいはさらに応用で使いやすくするために、本章でもう少しTight-bindingな近似を行い、解を具体化していきます。

# Tight-bindingな近似その2：重なり積分と飛び移り積分の近似

さて、改めて「強く束縛された状態」の描像に基づいた近似、つまり原子準位の広がりが小さいという近似を行っていきます。これにより簡単な場合であれば具体的なエネルギーバンド・固有関数が得られます。

### 重なり積分の近似

まずは**重なり積分**の近似を考えていきます。
具体的には、原子軌道$\varphi_m(\boldsymbol{r})$に対して、異なる格子点を中心に持つ次の積分を$0$と考えます。つまり

$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r} \simeq \delta_{m,n}\delta_{\boldsymbol{R},\boldsymbol{R}'}
$$

とします。この積分は前章で述べたようにOverlap Integral、重なり積分等と呼ばれるもので、前章では「小さい」という仮定の下進めていましたが、ここで思い切ってゼロとします。

この近似は、格子間距離に比べてWannier関数を近似している原子軌道の広がりが十分小さく、従って原子軌道同士の「重なり」が小さいという仮定に対応しています。

この近似により先ほど得た一般化固有値方程式において、行列$S_{\boldsymbol{k}}$が単位行列になり、

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

方程式は普通の固有値方程式

$$
\begin{Vmatrix}
(M_{\boldsymbol{k}})_{nm} - N\varepsilon_{n,\boldsymbol{k}}\delta_{nm}
\end{Vmatrix}
=0
$$

となり、問題は行列$M_{\boldsymbol{k}}/N$の対角化に帰着します。

### 飛び移り積分の近似

さらに行列$\left(M_{\boldsymbol{k}}\right)_{nm}$について、

$$
\begin{align*}
\left(M_{\boldsymbol{k}}\right)_{nm} 
&=
 \int \Phi_{n,\boldsymbol{k}}^*(\boldsymbol{r})\hat{H}\Phi_{m,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r},\\

&=
\sum_{\boldsymbol{R}'}
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}'-\boldsymbol{R})}\int
 \phi_n^*(\boldsymbol{r}-\boldsymbol{R}')
 \left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}''}V(\boldsymbol{r} - \boldsymbol{R}'')
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}\\
\end{align*}
$$



とできますが、ここで並進の対称性を利用して$N$個の$\boldsymbol{R}'$に対して積分は同じ値を取るので、$\boldsymbol{R}' = \boldsymbol{0}$の場合の$N$倍を考えることで

$$
\begin{align*}
&=
N
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R})}\int
 \phi_n^*(\boldsymbol{r})
 \left(
 -\frac{\hbar^2}{2m}\nabla^2
    +
   \sum_{\boldsymbol{R}''}V(\boldsymbol{r} - \boldsymbol{R}'')
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}\\

\end{align*}
$$

となります。さらに原子軌道関数が孤立原子のハミルトニアンの固有関数であることを利用し、また$\boldsymbol{R}'' = \boldsymbol{R}$とそれ以外の場合で分けることで

$$
\begin{align*}

 &=
N
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}
\varepsilon_{n}^{\rm a}
\int
 \phi_n^*(\boldsymbol{r})
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}

\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}''\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}'')
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}\\

 &=
N
\varepsilon_{n}^{\rm a}
\delta_{n,m}\\

&\>\>\>\>+
N

\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}''\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}'')
   \right)
 \phi_m(\boldsymbol{r})d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R})
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}\\

&\>\>\>\>+
N
\sum_{\boldsymbol{R}}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}}
\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}''\neq \boldsymbol{0}, \boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}'')
   \right)
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}\\
\end{align*}
$$

となります。

ここで、第2項の積分を

$$
\begin{align*}
\int
 \phi_n^*(\boldsymbol{r})
 \left(
   \sum_{\boldsymbol{R}''\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}'')
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
   \sum_{\boldsymbol{R}''\neq \boldsymbol{0}}V(\boldsymbol{r} - \boldsymbol{R}'')
   \right)
 \phi_m(\boldsymbol{r})d\boldsymbol{r}
 
 &\simeq
\Delta\varepsilon\delta_{nm}
 \end{align*}
$$

と定数の対角行列と近似されることもあるようです。このようにすればこの項はエネルギーを定数シフトさせるだけの効果となり、後の計算では無視していけるようになります。（本章では一応残しておきます。）

第3項は、原点と異なる格子点と、それと等しい局所ポテンシャルの積で、関数の局在性から近い位置の間の積分のみ値を持つと考え、特に原点から **隣接格子点**（多くの場合は第二隣接くらいまでらしい）へのベクトルを$\boldsymbol{R}_I$として

$$
\sum_{\boldsymbol{R}}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R})
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}

 \simeq
\sum_{\boldsymbol{R}_I}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I}
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_I)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d\boldsymbol{r}

$$

と近似されます。ここで$\sum_{\boldsymbol{R}_I}$は、隣接格子点のみ和を取ることを表します。

積分

$$
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_I)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d\boldsymbol{r}
 \equiv
 - t_{\boldsymbol{R}_I}^{n,m}
$$

は、上記のように符号を逆転させた$-t$で表示されることが多く（理由はよく知りませんが、同じ原子軌道間の積分を考えて原子軌道が原点に対して対称なら、（多分）重なりの大きい場所は同符号になり、それと負の値を持つ引力ポテンシャルの積の積分全体が負の値を持ちがちだからでしょうか）、また **"Transfer Integral"、"Hopping Integral"、「飛び移り積分」** 等と呼ばれます。（多分"t"は"transfer"の"t"だと思います）
多くの場合再隣接格子、時には次近接辺りまでの値が近似で採用されます。


また第4項は、全て中心が異なる関数の積分であり、同様にそれぞれの関数が局在していることから無視することができると考えます。（この形の積分を「3中心積分」等と呼ぶようです）

以上をまとめると、対角化すべき行列$M_{\boldsymbol{k}}/N$は、

$$
\begin{align*}
\left(M_{\boldsymbol{k}}\right)_{nm}/N
&\simeq

\varepsilon_{n}^{\rm a}
\delta_{n,m} - \Delta\varepsilon_{nm} - \sum_{\boldsymbol{R}_I}
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I}
 t_{\boldsymbol{R}_I}^{n,m}
\end{align*}

$$


となります。かなり簡略化されましたね。
さらに結晶場積分の値が全て等しいとして$\Delta\varepsilon\delta_{nm}$とおき、かつ展開に用いた軌道関数はエネルギーが近いものでしたので、思い切って原子準位$\varepsilon_{m_i}^{\rm a}$も全て等しいと近似し、その値を$\varepsilon^{\rm a}$と置けば、

$$
\begin{align*}
\left(M_{\boldsymbol{k}}\right)_{nm}/N
&\simeq

(\varepsilon^{\rm a} - \Delta\varepsilon) \delta_{n,m} -  \sum_{\boldsymbol{R}_I} 
e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I}
t_{\boldsymbol{R}_I}^{n,m}
\end{align*}

$$

となり、



単位行列$I$と、飛び移り積分の行列$T$から

$$
\frac{M_{\boldsymbol{k}}}{N} = (\varepsilon^{\rm a} - \Delta\varepsilon)I + 

T_{\boldsymbol{k}},\\

(T_{\boldsymbol{k}})_{nm} = - \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{n,m}
$$

で、最終的に固有値方程式

$$
\begin{Vmatrix}
(T_{\boldsymbol{k}})_{nm} - \varepsilon\delta_{nm}
\end{Vmatrix} = 0,\\
\varepsilon = \varepsilon_{\boldsymbol{k}} - (\varepsilon^{\rm a} - \Delta\varepsilon)
$$

となり、エネルギーバンドの形状は飛び移り積分$t_{\boldsymbol{R}_I}^{nm}$の値（と結晶構造）のみで決定されることになります。局所ポテンシャルと原子軌道関数はすでに得られているという仮定なので、後は飛び移り積分を計算して値を代入して対角化すれば解が得られますね。



## 原子軌道や分子軌道の名前を用いたバンドの呼び方

なお、この時、例えば展開に$s$軌道を選んだ場合、そこから得られるバンドを「$s$バンド」や「$s$的なバンド」、エネルギーが縮退した3つの$p$軌道で展開した場合に得られるバンドを「$p$バンド」「$p$的なバンド」などと呼びます。
また、（私はあまり詳しくないですが）展開の種類によっては分子軌道の用語「$\sigma$バンド」や「$\pi$バンド」なんかもよく見かけます。ここで「$\sigma$バンド」は$s, p_x, p_y$の展開、「$\pi$バンド」は残りの$p_z$軌道による展開に対応しているようです。


# Tight-bindingモデル

最後に、先ほど得られた固有値方程式

$$
\begin{Vmatrix}
(T_{\boldsymbol{k}})_{nm} - \varepsilon\delta_{nm}
\end{Vmatrix} = 0,\\
\varepsilon = \varepsilon_{\boldsymbol{k}} - (\varepsilon^{\rm a} - \Delta\varepsilon)
$$



を違う見方で見てみようと思います。

この方程式は飛び移り積分$t_{\boldsymbol{R}_I}$のみをパラメータとして持つ、エネルギーバンドを決定する方程式になっています。

そして、ここからは統計力学の範疇になるのですが、エネルギーバンドの形状（関数形$\varepsilon(\boldsymbol{k})$）が分かれば、様々な固体の物理量を計算することができるようになります。従ってエネルギーバンドは、固体の性質を端的に表現する関数であり、固体を表す「モデル」となります。

そこで今度は、唯一のパラメータ$t_{\boldsymbol{R}_I}$を計算で求めるのではなく、逆に$t_{\boldsymbol{R}_I}$を色々と変えたときに固体がどのような性質を持つのか、予言（予測）ができるようになります。
というわけで上記のような$t_{\boldsymbol{R}_I}$をパラメータとして表したエネルギーバンドを指して、Tight-binding Modelということも多いです。
（前章の出発点のように、固体のポテンシャルを1体ポテンシャルで近似＋固有状態を原子軌道で近似するのも、モデル化なので、そのような系を指してTight-bindingモデルということもあると思います）（この辺はだいぶ雰囲気で書いてるので全然違ってるかもしれません）

## 簡単な例でのエネルギーバンド＝Tight-bindingモデル

最後に、最も簡単な例として、1次元格子において、固有関数を$s$軌道を用いて展開した場合について具体的に考えてみます。$s$軌道は縮退がないので、もっとも粗い近似では1つの原子軌道（Bloch和）のみで、Wannier関数（Bloch関数）を展開できます。
この固有関数を$\varphi_{s,\boldsymbol{k}}(\boldsymbol{r})$と書きます。

なお、1次元系をどのように表現すればよいかちょっとややこしいのですが、最も手っ取り早い方法は基本格子ベクトル$\boldsymbol{a}_i$を、

$$
|\boldsymbol{a}_2| = |\boldsymbol{a}_3| = A \gg |\boldsymbol{a}_1|
$$

と置いてしまうことかと思います。これは実際の1次元系について、それ以外の方向との格子間距離が離れている状況を表しています。（2方向に結晶サイズ1、つまり$N_2 = N_3 =1$の周期的境界条件を課しても上手くいきそうな気もするのですが、なんかややこしいのでやめておきます。）こうすると後述するように、考慮しない方向の飛び移り積分が$0$となり、1次元的なバンドが得られます。

さて、格子間隔が開いたとはいえ、固有関数は先述のように一つの$s$軌道関数を用いて、

$$
\varphi_{s,\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\phi_s(\boldsymbol{r}-\boldsymbol{R})
$$

と表せます。ここでややこしいですが、格子は1次元でも、その中の（固有）関数は3次元的に広がっており、展開に用いる原子軌道も3次元関数です。

1つの軌道で展開した場合は行列$M_{\boldsymbol{k}}$は行列ではなくなり、固有値方程式は


$$
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} + \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}
 - \varepsilon_{s,\boldsymbol{k}}
=0
\\

\Rightarrow
\varepsilon_{s,\boldsymbol{k}}
=
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} + \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}
$$

と、ただの（？）方程式になります。
ここで$\varepsilon_{s}^{\rm a}$は$s$軌道の原子準位、$\Delta\varepsilon_{s}$は$s$軌道の結晶場積分、そして$\varepsilon_{s,\boldsymbol{k}}$が、電子の固有状態が持つ固有エネルギーの式（エネルギーバンド）を表します。


ここでここまではまだ3次元っぽいですが、ここで飛び移り積分$t_{\boldsymbol{R}_I}$について考えると、1次元系っぽくなります。
すなわち、$\boldsymbol{a}_2, \boldsymbol{a}_3$方向の再隣接格子へのベクトルを$\boldsymbol{A}$とし、その時の飛び移り積分を$-t_A$と置くと、

$$
t_A = - \int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} \pm \boldsymbol{A})\phi_s(\boldsymbol{r}\pm \boldsymbol{A})d\boldsymbol{r}
$$

となりますが、ここで$A$は大きい数と置いたことから、$\boldsymbol{a}_2, \boldsymbol{a}_3$方向の飛び移り積分は$t_A\simeq0$となります。この結果エネルギーバンドは

$$
\varepsilon_{s,\boldsymbol{k}}
\simeq
\varepsilon_{s}^{\rm a} - \Delta\varepsilon_{s} - \sum_{R_1 = \pm a}e^{ik_1 R_1} t_{\boldsymbol{R}_I}
$$

と、波数空間における1次元関数（$k_2, k_3$方向には定数関数）のエネルギーバンドが得られます。

さらに、$V(\boldsymbol{r})$が球対称ポテンシャルを仮定しており、また$s$軌道も角度部分は定数、つまり球対称であるため、$R_I = \pm a$いずれも同じ積分の値を持ちます。

また、その値を以下のように$-t$と置くと（今$s$軌道関数は球対称関数のため、$\boldsymbol{a}$離れた関数の間の重なりが大きい領域は同符号であると考えて良さそうです。従って飛び移り積分の値は負（$t>0$）になります。）

$$
-\int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} + \boldsymbol{a}_1)\phi_s(\boldsymbol{r} + \boldsymbol{a}_1)
=

-\int \phi_s^*(\boldsymbol{r})V(\boldsymbol{r} - \boldsymbol{a}_1)\phi_s(\boldsymbol{r} - \boldsymbol{a}_1)
\equiv t >0
$$


1次元関数の固有エネルギー（エネルギーバンド）は

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

飛び移り積分はざっくり、格子間距離に対する、孤立原子ポテンシャルの強さ（広がり）と、原子軌道の広がりに依存するのでした。従って、上記のTight-bindingモデルは、原子間距離が近かったり、あるいは孤立原子ポテンシャルが強ければ、バンドの幅が大きくなり、逆に格子間距離が広かったり、孤立原子ポテンシャルが弱ければ、バンドの幅が狭くなることを表しています。^[極端な場合、格子間距離が大きく離れていて、飛び移り積分がゼロの場合はどうなるでしょうか？その場合はバンド幅がゼロ、つまりBZ内のBloch波数$\boldsymbol{k}$の数の、縮退した固有状態を持ってる状態です。これは、結局N個の原子が孤立していて、そのどこか一つに電子が束縛されているような状態を意味しています。]で、バンドの幅が広かったり、狭かったりするとどうなるのか、という辺りは固体物理の教科書に譲ることにします。

# おわりに

本章ではTight-binding近似、そしてTight-bindingモデルの導出の後編として、
重なり積分が小さく、


$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r} \simeq \delta_{m,n}\delta_{\boldsymbol{R},\boldsymbol{R}'}
$$

と近似し、さらに隣接格子間の飛び移り積分

$$
\int
 \phi_n^*(\boldsymbol{r})
 
   V(\boldsymbol{r} - \boldsymbol{R}_I)
 
 \phi_m(\boldsymbol{r}-\boldsymbol{R}_I)d\boldsymbol{r}
 \equiv
 - t_{\boldsymbol{R}_I}^{n,m}
$$
のみを考えたり、その他色々して、最終的にエネルギーバンドを求める固有値方程式

$$
\begin{Vmatrix}
(T_{\boldsymbol{k}})_{nm} - \varepsilon\delta_{nm}
\end{Vmatrix} = 0,\\
\varepsilon = \varepsilon_{\boldsymbol{k}} - (\varepsilon^{\rm a} - \Delta\varepsilon),\\
(T_{\boldsymbol{k}})_{nm} = - \sum_{\boldsymbol{R}_I}e^{i\boldsymbol{k}\cdot\boldsymbol{R}_I} t_{\boldsymbol{R}_I}^{n,m}
$$

を導きました。

本章を通して、Tight-binding近似のもと、

- 固有関数をどのように表すことができるか
- 固有関数が従うべき方程式をどのように近似（簡略化）していけるか

が、分かったかと思います。

いよいよ残すは上記の理解をベースに、第二量子化表示のTight-bindingモデルを導けば本稿の目標は達成です。

・・・が、その前に次章では、前章と本章で出てきた「飛び移り積分」「重なり積分」の物理的な意味を少し~~妄想して~~考えてみようと思います。（なぜなら学生時代の私がずっと疑問に思っていたからです）
