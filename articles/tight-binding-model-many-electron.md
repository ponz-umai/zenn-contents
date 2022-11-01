---
title: "多電子状態の波動関数"
emoji: "🌟"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# はじめに
前回の[水素原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-hydrogen-atom)では「1つの原子核＋1つの電子」の場合について扱い、電子の固有状態を求めました。

最終目標は「多数の原子核（や分子）＋多数の電子」からなる固体で、以下のような状態について理解することです。（抽象画みたいになってしまいましたが・・・）

![](/images/tb/many-ele-atm.png)

多数の原子核の相互作用を受けながら、多数の電子がそれぞれに相互作用をしあっている状態です。わけわからんですね。ただ一つ一つ理解していけば大丈夫です。（図はいずれ、もっときれいにします。。）

「多数の原子核＋多数の電子」からなる固体の電子状態に近づくためには原子核を増やすか、電子を増やすかなのですが、話の流れ上まずは電子を増やして「1つの原子核＋多数の電子」の方向に進むことにします。

描像としては以下のイメージです。

![](/images/tb/many-electron.png)

水素原子の場合と同様に原子核を原点に選び静止していると仮定します。原子核に束縛された複数の電子を考え、それぞれの電子の座標を$\boldsymbol{r} _i$、$i$番目の電子から$j$番目の電子の位置ベクトルを$\boldsymbol{r} _{ij}$と置きます。

原子核と電子$i$の間の（引力）クーロン相互作用ポテンシャルを$V_\mathrm{a}(r_i)$（atomのa）, 電子間の（斥力）クーロン相互作用ポテンシャルを$V_\mathrm{e}(r_{ij})$とします。

この時シュレディンガー方程式は

$$
\mathcal{H}\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) = E\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots),\\
\mathcal{H} = \sum_i\left(  -\frac{\hbar^2}{2m}\nabla_i^2\right) + \sum_iV_a(r_i) + \frac{1}{2}\sum_{i\neq j}V_{\mathrm{e}}(r_{ij})
$$

となります。ここで$\nabla_i$は$i$番目の粒子に作用する微分演算子を表します。また$\mathcal{H}$の第3項では二重に数えた分を戻すために2で割っています。

早速この方程式を解いていきたいのですが、その前段としてまず多粒子系を扱う際の特別な性質を理解する必要があり、そちらが結構長くなってしまったので「1原子核に束縛された多電子」の問題は次章に持ち越すことにして、本章では多電子系の性質についてまとめることとします。

# 多粒子系の量子力学

さて、上記シュレディンガー方程式では水素原子のシュレディンガー方程式の場合と比べて
- ハミルトニアンを筆記体$\mathcal{H}$で書く
- 波動関数を大文字$\Phi$で書き、座標$\boldsymbol{r}_i$の多変数関数とする
- エネルギーも多粒子のエネルギーの合計という気持ちを込めて大文字$E$で表す

と、記号を変更しました。というのも、粒子数が複数の場合には特殊な性質が色々と出現するので一粒子の場合と区別する必要があるからです。また後ほど出てきますが、多粒子系も最終的には一粒子系の方程式に帰着するので、その移り変わりを明確にすることで「今どっちを考えているのか」がわかりやすくなる効果もあります。

多粒子系の取り扱いについても多数の素晴らしい書籍やWEB上の講義ノートが存在しますので、あえて私がその劣化版を書き連ねることもないかと思いつつも、自分の復習や将来読み返すために、必要な内容を最小限に絞りつつまとめていこうと思います。

## 多粒子のシュレディンガー方程式

多数の粒子の状態を取り扱うにあたり出てくる特殊事情について、まずは波動関数ですが、多粒子系のハミルトニアン演算子が

$$\mathcal{H} = \sum_i\left(  -\frac{\hbar^2}{2m}\nabla_i^2\right) + \sum_iV_a(r_i) + \frac{1}{2}\sum_i\sum_{j\neq i}V_{\mathrm{e}}(r_{ij})
$$

と複数の位置座標に作用する演算子であることから、その解となる波動関数も

$$
\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots)
$$

と、多数の位置座標を引数として持つ多変数関数となります。これは
1番目の粒子が座標$\boldsymbol{r}_1$、2番目の粒子が座標$\boldsymbol{r}_2$、・・・で観測される確率（密度）」というように全ての粒子の状態をひとまとめにした関数です。

とは言っても、1粒子のハミルトニアンも3次元座標$(x,y,z)$や$(r,\theta,\phi)$に対する演算子で、その解である波動関数も$\varphi(x,y,z)$と多変数関数であったので、その意味では事情は同じです。

そして一粒子の場合（例えば水素原子の場合）ハミルトニアンを$r,\theta,\phi$に対する演算子に分解し、変数分離法を使って波動関数も$R(r),\Theta(\theta),\Phi(\phi)$と分解したように、多くの場合様々な近似法を駆使して、多粒子のハミルトニアン・波動関数も1粒子演算子・波動関数に分解して解くことになります。

まず波動関数の性質を理解していくことを目標に、一体演算子への近似法は次章以降で扱うこととして、本章ではシュレディンガー方程式が変数分離できる形をまず想定して解の性質を理解していくこととします。

:::message alert
本章ではスピンについては扱わず、次章でスピンを考慮した取り扱いについて記載する予定です。
:::

## 一体演算子の和からなるシュレディンガー方程式

つまりハミルトニアンが

$$
\mathcal{H} \simeq \sum_i\hat{H}_i
$$

と分解できた場合について取り扱います。これは粒子間の相互作用がない状態を意味します。実際は、粒子間の相互作用をうまく近似して、1粒子のみが関わるポテンシャルに書き換えた状態に対応します。

このときシュレディンガー方程式は

$$
\sum_i\hat{H}_i\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) = E \Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots)
$$

となり、この微分方程式は

$$
\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) =\phi_1(\boldsymbol{r}_1)\phi_2(\boldsymbol{r}_2)\cdots\phi_i(\boldsymbol{r}_i)\cdots
$$

と置いて両辺を$\phi_1(\boldsymbol{r}_1)\phi_2(\boldsymbol{r}_2)\cdots\phi_i(\boldsymbol{r}_i)\cdots$で割ることで、

$$
\frac{\hat{H}_1\phi_1(\boldsymbol{r}_1)}{\phi_1(\boldsymbol{r}_1)} + \frac{\hat{H}_2\phi_2(\boldsymbol{r}_2)}{\phi_2(\boldsymbol{r}_2)} + \cdots = E
$$

となり、左辺の各項を

$$
\frac{\hat{H}_i\phi_i(\boldsymbol{r}_i)}{\phi_i(\boldsymbol{r}_i)}  = \epsilon_i\\
\Rightarrow\hat{H}_i\phi_i(\boldsymbol{r}_i)  = \epsilon_i\phi_i(\boldsymbol{r}_i)
$$

と置くことで最終的には一粒子のシュレディンガー方程式を解き、その固有関数を求めればよいことになります。

## 多粒子系の波動関数の性質

### 反対称性
上記の変数分離により、シュレディンガー方程式

$$
\sum_i\hat{H}_i\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) = E \Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots)
$$

の解は、それぞれのシュレディンガー方程式を解いて得られた固有関数を
$\varphi_\lambda, \varphi_\mu,\cdots,\varphi_\xi,\cdots$として

$$
\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) =\varphi_\lambda(\boldsymbol{r}_1)\varphi_\mu(\boldsymbol{r}_2)\cdots\varphi_\xi(\boldsymbol{r}_i)\cdots
$$

となります。数学的にはこれでOK！なのですが、ここで新たな条件が付け加わることになります。

:::message
複数の同種粒子を互いに区別することはできない
:::

これは「同種粒子の不可弁別性」と呼ばれる性質で、平たく言うと同じ種類の粒子（例えば電子）が複数個あるとき、ある電子を観測した際にその電子が「電子1」なのか、「電子2」なのか・・・を判別することができない。という性質です。すなわち私たちがわかるのは、「座標$\boldsymbol{r}_i$に電子が観測された」という事実だけ、ということです。実際、電子に印をつけることはできませんし、古典力学の場合のように電子の移動を逐次追っていくことはできませんので、観測された電子がどの電子なのかを特定することはできない、というのは納得感はあります。

とはいえだからどうした？という感じで、実際電子はそれぞれ個別に存在しているのではないのか？と思いますが、この条件を正しく取り入れることで様々な物理的な描像を正しく理解することができるようになるようです。つまり上述の条件は便宜的なものではなく、何かミクロの世界の本質を表しているもののようです。（この辺りはまだ全然良くわからないので、もし宝くじに当たったらちゃんと研究してみたいなーなどという願望を抱いています。）

さて、上記の条件を、電子の状態を表す波動関数に上手く取り込む必要があります。それには次のように考えます。

:::message
多粒子の分布を表す波動関数$\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots)$は、その座標変数を入れ替えても結果が変わらない（座標変数の置換対称性を持つ）
:::

今$\boldsymbol{r}_i(\boldsymbol{r}_j)$は、「粒子$i(j)$の座標」と定義しており、

$$
\Phi(\cdots\boldsymbol{r}_i, \boldsymbol{r}_j,\cdots)
$$

は、「$i$番目の粒子が座標$\boldsymbol{r}_i$、$j$番目の粒子が座標$\boldsymbol{r}_j$で観測される確率」を表しています。それを入れ替えても結果が変わらないということは、「$i$番目の粒子が座標$\boldsymbol{r}_i$、$j$番目の粒子が座標$\boldsymbol{r}_j$にいる（観測される）状態」と「$i$番目の粒子が座標$\boldsymbol{r}_j$、$j$番目の粒子が座標$\boldsymbol{r}_i$にいる（観測される）状態」が等しい、ということです。分かったような分らんような気分ですが、一応粒子の不可弁別性を上手く表せているように思います。

実際には粒子の存在確率は波動関数の絶対値の二乗なので、

$$
\left|\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_i,\cdots)\right|^2 
= \left|\Phi(\boldsymbol{r}_2, \boldsymbol{r}_1,\cdots,\boldsymbol{r}_i,\cdots)\right|^2\\
=\cdots=\left|\Phi(\boldsymbol{r}_i, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_1,\cdots)\right|^2=\cdots
$$

を満たす、という条件が新たに課されることになります。絶対値の二乗を外すと、絶対値の二乗で1になる係数$C, |C|^2=1$を用いて

$$
\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_i,\cdots,\boldsymbol{r}_j,\cdots)=
C\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_j,\cdots,\boldsymbol{r}_i,\cdots)
$$

と書けます。ここで電子は、

$$
C = -1
$$

であるようです。すなわち「波動関数の座標を入れ替えるとマイナス1がかかる」という性質を持ちます。このような粒子のことを「フェルミ粒子」または「フェルミオン」と呼び、その性質を波動関数の「反対称性」と呼びます。なお、座標の入れ替えに対する係数$C＝1$であるような粒子を「ボース粒子」「ボソン」と呼びます（ボゾンとも）。

本稿ではおそらく扱いませんが、例えばヘリウム原子や（本稿ではずっと無視されることになる）原子核の振動を量子化したフォノン、光子（フォトン）等がボース粒子の性質を持ち、特に顕著な効果として液体ヘリウムの超流動（ボース・アインシュタイン凝縮）等へつながっていたりします。

~~Cは一般に複素数でいいはずですが、この係数が複素数になるような粒子はあるのでしょうか？という疑問を抱いてしまいますが一旦忘れて進みます。~~


## スレーター行列式

さて話を戻して、シュレディンガー方程式$\sum_i\hat{H}_i\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) = E \Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots)$は普通に変数分離で解くと$\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) =\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) =\varphi_\lambda(\boldsymbol{r}_1)\varphi_\mu(\boldsymbol{r}_2)\cdots\varphi_\xi(\boldsymbol{r}_i)\cdots$となりました。しかしこの解は反対称性を満たしません。どうすればよいでしょうか。

ここで線形微分方程式の性質を利用します。すなわち、

:::message
線形微分方程式の解の線形結合もまた、元の微分方程式の解となる
:::

という性質を利用し、変数分離解$\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots) =\varphi_\lambda(\boldsymbol{r}_1)\varphi_\mu(\boldsymbol{r}_2)\cdots\varphi_\xi(\boldsymbol{r}_i)\cdots$の線形結合を上手くとることで反対称性を満たす解を作ります。

いきなり一般の多粒子状態を扱うとややこしくなるので、例として2粒子状態のシュレディンガー方程式

$$
\left(\hat{H}_1+\hat{H}_2\right)\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2) = E \Phi(\boldsymbol{r}_1, \boldsymbol{r}_2)
$$

を解いて変数分離解

$$
\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2) =\varphi_a(\boldsymbol{r}_1)\varphi_b(\boldsymbol{r}_2)
$$

を得たとします。この時波動関数を

$$
\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2) =\frac{1}{\sqrt{2}}\left\{
    \varphi_a(\boldsymbol{r}_1)\varphi_b(\boldsymbol{r}_2) - \varphi_b(\boldsymbol{r}_1)\varphi_a(\boldsymbol{r}_2)\right\}
$$

と書くと（$1/\sqrt{2}$は規格化係数）、確かに座標の入れ替えに対して

$$
\begin{align*}
\Phi_\mathrm{F}(\boldsymbol{r}_2, \boldsymbol{r}_1) 
&=\frac{1}{\sqrt{2}}\left\{
    \varphi_a(\boldsymbol{r}_2)\varphi_b(\boldsymbol{r}_1) - \varphi_b(\boldsymbol{r}_2)\varphi_a(\boldsymbol{r}_1)
    \right\}\\
&=\frac{1}{\sqrt{2}}\left\{
    -\varphi_a(\boldsymbol{r}_1)\varphi_b(\boldsymbol{r}_2) +\varphi_a(\boldsymbol{r}_1)\varphi_a(\boldsymbol{r}_2)
    \right\}\\
&=-\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2)
\end{align*}
$$

と、反対称になっていることがわかります。
少し先取りして書くと反対称化した波動関数$\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2) =\frac{1}{\sqrt{2}}\left\{    \varphi_a(\boldsymbol{r}_1)\varphi_b(\boldsymbol{r}_2) - \varphi_b(\boldsymbol{r}_1)\varphi_a(\boldsymbol{r}_2)\right\}$は行列式を用いて

$$
\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2) =\frac{1}{\sqrt{2}}
\begin{vmatrix}
\varphi_a(\boldsymbol{r}_1) & \varphi_b(\boldsymbol{r}_1) \\
\varphi_a(\boldsymbol{r}_2) & \varphi_b(\boldsymbol{r}_2) \\
\end{vmatrix}
$$

と書くことができます。N粒子の反対称化した波動関数も同様に、

$$
\Phi_\mathrm{F}(\cdots,\boldsymbol{r}_i, \boldsymbol{r}_j,\cdots,\boldsymbol{r}_N) = 
-\Phi_\mathrm{F}(\cdots,\boldsymbol{r}_j, \boldsymbol{r}_i,\cdots,\boldsymbol{r}_N) =\cdots
$$

と座標の入れ替えごとにマイナスが掛かるような線形結合で表したいのですが、これも上述のように「行列の行を入れ替えるとマイナスがかかる」という行列式の性質を用いることで、
変数分離で得られた固有関数を$\varphi_\lambda, \varphi_\mu,\cdots,\varphi_\xi$としてそれらを並べた行列式

$$
\begin{align*}
\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N) &= 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_\lambda(\boldsymbol{r}_1) & \varphi_\mu(\boldsymbol{r}_1) & \cdots & \varphi_\xi(\boldsymbol{r}_1)\\
\varphi_\lambda(\boldsymbol{r}_2) & \varphi_\mu(\boldsymbol{r}_2) & \cdots & \varphi_\xi(\boldsymbol{r}_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\boldsymbol{r}_N) & \varphi_\mu(\boldsymbol{r}_N) & \cdots & \varphi_\xi(\boldsymbol{r}_N)
\end{vmatrix}\\
&=\frac{1}{\sqrt{N!}}\sum_P(-1)^P\hat{P}\varphi_\lambda(\boldsymbol{r}_1)\varphi_\mu(\boldsymbol{r}_2)\cdots\varphi_\xi(\boldsymbol{r}_i)\cdots
\end{align*}
$$

で書くことができます。この時座標$\boldsymbol{r}_i$と座標$\boldsymbol{r}_j$の入れ替えは、行$i$と行$j$の入れ替えに対応し、行列式の性質から1回の入れ替えごとに波動関数の符号が変わる、という性質を表現できています。

また第二式は行列式の別の表現で、$\hat{P}$はその右側にある関数の間の座標を置換する演算子で、$(-1)^P$は置換1回ごとに$(-1)$が掛かる（偶置換なら正、奇置換なら負）という意味です。そして$\sum_P$は全ての置換について和を取ることを表します。

2粒子系で書き下すと、

$$
\sum_P(-1)^P\hat{P}\varphi_a(\boldsymbol{r}_1)\varphi_b(\boldsymbol{r}_2) = 
    \varphi_a(\boldsymbol{r}_1)\varphi_b(\boldsymbol{r}_2) - \varphi_a(\boldsymbol{r}_2)\varphi_b(\boldsymbol{r}_1)
$$

第二項が$a$と$b$の順番が違いますが、先ほどの式と同じですね。第一項はゼロ置換なので$(-1)^P = 1$, 第二項は一回置換しているので$(-1)^P = -1$です。

### パウリの排他原理

このようにしてようやく、多体の波動関数を得ることができました。

ここで波動関数をスレーター行列式で表した際、その解は以下の重要な性質を持ちます。

:::message 
スレーター行列式に含まれる1粒子固有関数$\varphi_\lambda, \varphi_\mu,\cdots,\varphi_\xi,\cdots$はすべて異なる関数でなければならない
:::

これは「パウリの排他原理」と呼ばれる性質で、行列式の性質から直接導かれます。すなわち行列式の同じ行に二つ以上の同じ関数があれば、行列の列を入れ替えると行の場合と同様に符号が入れ替わるという性質を用いて、2列目と最終列を入れ替えます。（見た目は同じですが入れ替えたことで負号がかかっています）

$$
\begin{align*}
\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N) &= 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_\lambda(\boldsymbol{r}_1) & \varphi_\mu(\boldsymbol{r}_1) & \cdots & \varphi_\mu(\boldsymbol{r}_1)\\
\varphi_\lambda(\boldsymbol{r}_2) & \varphi_\mu(\boldsymbol{r}_2) & \cdots & \varphi_\mu(\boldsymbol{r}_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\boldsymbol{r}_N) & \varphi_\mu(\boldsymbol{r}_N) & \cdots & \varphi_\mu(\boldsymbol{r}_N)
\end{vmatrix}\\
&=-\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_\lambda(\boldsymbol{r}_1) & \varphi_\mu(\boldsymbol{r}_1) & \cdots & \varphi_\mu(\boldsymbol{r}_1)\\
\varphi_\lambda(\boldsymbol{r}_2) & \varphi_\mu(\boldsymbol{r}_2) & \cdots & \varphi_\mu(\boldsymbol{r}_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\boldsymbol{r}_N) & \varphi_\mu(\boldsymbol{r}_N) & \cdots & \varphi_\mu(\boldsymbol{r}_N)
\end{vmatrix}\\
&=-\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N)\\
&\Rightarrow \Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N) = 0
\end{align*}
$$

と示せます。

この性質はかなり重要で、多粒子系の状態を考える際は必ず考えることになります。


後述もしますが、化学の周期表を勉強した際に、ソロバンのように電子を軌道に一つづつ"詰めて”いくように習ったかと思います。その頃は疑いもせずに受け入れていた描像ですが、あれも「パウリの排他原理があるため、電子が一つ増えるごとに各電子はそれぞれ異なる固有状態を取っている」という理屈なのだとわかります。また、こちらも先取りですが固体内の電子の固有状態は「エネルギーバンド」というグネグネした道に下から電子を埋めていくように考えるのですが、それもこの原理から来ています。

と言っても抽象的な記号が続いたので（$\varphi_\xi$とか）、改めて多粒子系の性質をまとめた上で、具体的な問題設定を扱って本章を終わりにしようと思います。

## ここまでのまとめと具体例

というわけで、多粒子系を扱う場合でも粒子間の相互作用を扱わず、ハミルトニアンが一粒子の演算子の和で書けた場合は、波動関数は

$$
\mathcal{H} \simeq \sum_i\hat{H}_i
$$

で書けた場合は、それぞれ一体のシュレディンガー方程式

$$
\hat{H}_i\phi_i(\boldsymbol{r}_i)  = \epsilon_i\phi_i(\boldsymbol{r}_i)
$$

を解いて固有関数を求め、（フェルミオンの場合は）得られた固有関数の積を

$$
\begin{align*}
\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N) &= 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_\lambda(\boldsymbol{r}_1) & \varphi_\mu(\boldsymbol{r}_1) & \cdots & \varphi_\xi(\boldsymbol{r}_1)\\
\varphi_\lambda(\boldsymbol{r}_2) & \varphi_\mu(\boldsymbol{r}_2) & \cdots & \varphi_\xi(\boldsymbol{r}_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\boldsymbol{r}_N) & \varphi_\mu(\boldsymbol{r}_N) & \cdots & \varphi_\xi(\boldsymbol{r}_N)
\end{vmatrix}\\
&=\frac{1}{\sqrt{N!}}\sum_P(-1)^P\hat{P}\varphi_\lambda(\boldsymbol{r}_1)\varphi_\mu(\boldsymbol{r}_2)\cdots\varphi_\xi(\boldsymbol{r}_i)\cdots
\end{align*}
$$

と反対称化することで多粒子系の波動関数が得られることがわかりました。

また重要な性質として、「スレーター行列式に含まれる1粒子固有関数$\varphi_\lambda, \varphi_\mu,\cdots,\varphi_\xi,\cdots$はすべて異なる関数でなければならない」というパウリの排他原理を持つことがわかりました。

上記を踏まえて、具体的にN個の自由電子のハミルトニアン

$$
\mathcal{H} = \sum_{i=1}^N\left( \frac{-\hbar^2}{2m}\nabla_i^2 \right)
$$

のみを考えた場合のシュレディンガー方程式を解いてみます。
上でまとめたように、上記ハミルトニアンは一体のハミルトニアンに分解でき、かつ全て同じ形になるので、ただ一つのシュレディンガー方程式

$$
\frac{-\hbar^2}{2m}\nabla^2 \phi(\boldsymbol{r}) = \epsilon\phi(\boldsymbol{r})
$$

を解けばよく、一体のシュレディンガー方程式の固有関数と固有エネルギーは、ラベル$\boldsymbol{k}$を用いて

$$
\varphi_{\boldsymbol{k}} = \frac{1}{\sqrt{V}}e^{i\boldsymbol{k}\cdot\boldsymbol{r}},\\
\epsilon_{\boldsymbol{k}} = C,\\
\boldsymbol{k} = (k_x, k_y, k_z), k_i = \frac{2\pi}{L}n_i,\\
(i = x,y,z, n_i = 0,\pm1,\pm2,\cdots)
$$

と書けます。ここで境界条件として一辺$L$の立方体の周期的境界条件を考えました。

そこで多粒子の波動関数は、上記の固有関数のうち**すべて異なる関数**をN個選んで作ったスレーター行列式

$$
\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N) = 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_{\boldsymbol{k}_a}(\boldsymbol{r}_1) & \varphi_{\boldsymbol{k}_b}(\boldsymbol{r}_1) & \cdots & \varphi_{\boldsymbol{k}_n}(\boldsymbol{r}_1)\\
\varphi_{\boldsymbol{k}_a}(\boldsymbol{r}_2) & \varphi_{\boldsymbol{k}_b}(\boldsymbol{r}_2) & \cdots & \varphi_{\boldsymbol{k}_n}(\boldsymbol{r}_2)\\
& \cdots & \cdots\\
\varphi_{\boldsymbol{k}_a}(\boldsymbol{r}_N) & \varphi_{\boldsymbol{k}_b}(\boldsymbol{r}_N) & \cdots & \varphi_{\boldsymbol{k}_n}(\boldsymbol{r}_N)
\end{vmatrix}
$$

で表されます。特に最低エネルギー状態（基底状態）の波動関数は、**固有エネルギーが小さい順からN個**選んだスレーター行列

$$
\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N) = 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_{\boldsymbol{k}_0}(\boldsymbol{r}_1) & \varphi_{\boldsymbol{k}_1}(\boldsymbol{r}_1) & \cdots & \varphi_{\boldsymbol{k}_N}(\boldsymbol{r}_1)\\
\varphi_{\boldsymbol{k}_0}(\boldsymbol{r}_2) & \varphi_{\boldsymbol{k}_1}(\boldsymbol{r}_2) & \cdots & \varphi_{\boldsymbol{k}_N}(\boldsymbol{r}_2)\\
& \cdots & \cdots\\
\varphi_{\boldsymbol{k}_0}(\boldsymbol{r}_N) & \varphi_{\boldsymbol{k}_1}(\boldsymbol{r}_N) & \cdots & \varphi_{\boldsymbol{k}_N}(\boldsymbol{r}_N)
\end{vmatrix}
$$


で書かれます。かなり大雑把に書いていますが、具体的には一番エネルギーが小さいのは

$$
k_x = k_y = k_z = 0\Rightarrow\epsilon_{\boldsymbol{k}} = 0
$$

で、その次が

$$
k_x = \frac{2\pi}{L}\times 1,  k_y = k_z = 0,\\
\cdots,\\
\epsilon_{\boldsymbol{k}} = \frac{2\pi^2\hbar^2}{mL^2}
$$
の3通りの状態です。

固有エネルギーは一粒子シュレディンガー方程式を解いて得られた固有エネルギーを小さい順に足し合わせた

$$
E = \epsilon_{\boldsymbol{k}_0} + \epsilon_{\boldsymbol{k}_1} + \cdots + \epsilon_{\boldsymbol{k}_N}
$$

となります。

ここまで考えれば、多粒子系のシュレディンガー方程式・波動関数の取り扱いも親しみを感じてくるのではないでしょうか。

::: message alert
ここではスピンを考えていませんが、スピンを考えて固有状態をしたから"詰めて"行くことで固体物理の教科書の最初の方に乗っている、「フェルミ球」「フェルミエネルギー」の式と一致するようになります。

次回改めてスピンについて取り扱い、そこで再度この問題を取り上げようと思います。
:::

# おわりに

今回は1体のハミルトニアンのみを題材に、多粒子系の波動関数が満たす性質「反対称性」を扱い、その帰結としてパウリの排他原理が現れることを見ました。

この効果によって、自由粒子のシュレディンガー方程式で見たように、基底状態の多体の波動関数は「エネルギーの小さい順から一体の波動関数を"詰めて"行く」ようにして構成されることが分かりました。

ただし本章ではフェルミ粒子の性質「スピン」を扱っておらず、もう少し補足が必要な状態です。
次章ではミクロな粒子が持つ謎の性質「スピン」を導入し、今回取り扱った電子の性質に説明を付け加える予定です。

その先はいよいよ多電子原子を扱っていく予定です。そこでは先述のように（高校の化学などで習うように）原子番号が増える順番に、電子がひとつづつ詰まっていく描像が、量子力学では上記パウリの排他律を考慮した波動関数の構成に対応していることがわかります。ここを学んだ際は慣れ親しんだ原子のイメージが突然数式によって説明され、個人的にとても印象深い体験として記憶に残っています。

というわけで今回はこんなところで。もっと端的にテンポよく書いていきたいのですがついつい長くなってしまうのが悩みどころです。