---
title: "電子のスピンを考慮した多電子系の波動関数"
emoji: "🗂"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# はじめに
前章[多電子状態の波動関数](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron)では多電子状態の波動関数が満たすべき「反対称性」を導入し、その要件がスレーター行列式で満たされることを見ました。そしてその結果「スレーター行列式に含まれる1粒子固有関数はすべて異なる関数でなければならない」というパウリの排他原理が導かれました。

その結果基底状態の波動関数は、エネルギーが低い順に一電子固有状態を"詰めて"行った構造となりました。

ところがここで一点補足があり、電子の状態が従う関数としてスピンなるものも考える必要があります。

電子の状態を表す波動関数として、これまでは位置座標$\boldsymbol{r}$の関数$\varphi(\boldsymbol{r})$を考え、その意味は「位置座標$\boldsymbol{r}$で電子が観測される確率密度」（正確には絶対値の二乗を取ったもの）でした。

しかし実は同じ位置座標の関数で表される状態$\boldsymbol{r}$の関数$\varphi(\boldsymbol{r})$は、上向きスピン状態と下向きスピン状態の二つの状態を取り得ることが分かっています。

とはいえこれだけだとなんのこっちゃという感じなので、本章で電子の持つ新たな性質「スピン」についてまとめていきます。

なお名前の通り「スピン」は回転、すなわち角運動量に対応した性質を持っており、普通の量子力学の教科書では電子の（軌道）角運動量の話題の後に「ちなみに」みたいな感じで導入されることが多いですが、こと「Tight-bindingについて理解する」という目標に限れば（多分）角運動量について触れずに行けると思っているので、本稿ではあくまで「波動関数の新たなラベル（座標）」という位置づけに限ってスピンを導入していこうと思います。

余力がどこかで生まれれば、別ページで角運動量演算子から始めてちゃんと書こうと思います。（いつか追記予定）

# スピンの必要最低限の導入

## 実験事実

最大限に省略して書くと、電子の様子を注意深く観測したところ、電子は3次元空中の運動状態によらない謎の角運動量を持っていることがわかり、かつその値はある軸を中心に丁度同じ大きさで「正の角運動量」「負の角運動量」（不正確だけど「左回転」「右回転」みたいなもん）に対応する2つの値のみを持つらしいことが実験から分かったようです。

（今のところの私の知る限り「どうやらそうらしい」としか言えないのでこのような自信なさげな書き方になっております。）

## スピン変数とスピン関数
したがって、今までは一粒子の固有状態を、ラベル（例えば自由電子なら$\boldsymbol{k} = (k_x,k_y,k_z)$、水素原子なら$n,l,m$）で分類して、$\varphi_1, \varphi_2,\cdots$と書いていましたが、その固有状態がさらに「正の角運動量を持つ状態」「負の角運動量を持つ状態」の両方の状態を持ちうる、ということになります。

### ラベル的な表現

このような二つの状態を、ラベル$n$（水素原子中の電子のように$n,m,l$等と複数含む場合も含む）で表される座標の関数$\varphi_n(\boldsymbol{r})$（**軌道部分**とか**軌道関数**とか呼んだりします。以降この呼び方を使います）に、
「正の角運動量状態」を上向き矢印$\uparrow$や$+$、「負の角運動量状態」を下向き矢印$\downarrow$や$-$でラベルを付けて、

- 正の角運動量状態：$\varphi_{n\uparrow}, \varphi_{n+}$
- 負の角運動量状態：$\varphi_{n\downarrow},\varphi_{n-}$
- 一般の場合：$\varphi_{n\sigma}$

と書いたり、正の角運動量状態はそのままで、負の角運動量状態を$\overline{\varphi_n}$と書いたりします。

矢印$\uparrow, \downarrow$と対応して、「正の角運動量状態」を「上向きスピン状態」、「負の角運動量状態」を「下向きスピン状態」と呼びます。本稿でも以降、そのように書きます。

### スピン変数とスピン関数を用いた表現

大抵の場合は上記の表現で乗り切れるのですが、例えば「軌道部分が複数の状態の重ね合わせ状態で、かつスピン部分も「上向き」と「下向き」の重ね合わせの状態」みたいな関数を考えるには上記の表現では限界があります。
（実際、そのような例を次章で扱う予定です）

そこで一般的には、座標を変数とした波動関数$\varphi(\boldsymbol{r})$の解釈「位置$\boldsymbol{r}$で電子が観測される確率密度」（正確には絶対値の二乗を取ったもの）と対応させて、「スピン変数」$\sigma$を新たな波動関数の引数として加え、「位置$\boldsymbol{r}$スピン$\sigma$で電子が観測される確率密度」

$$
\varphi(\boldsymbol{r},\sigma)
$$

を電子の従う波動関数と考えることになります。
本稿ではその時々で便利な表現を使い分けることにします。

::: message alert
ここまでで、あるいはこれからの説明を読む中で、「スピンは状態なのか変数（座標）なのかラベル（量子数）なのかどっちなんだ？」と混乱するかもしれません。というか私もいまだに混乱しています。

さしあたっては、電子は座標や運動状態とは独立に、謎の2種類の角運動量を持った状態を持つらしい、ということを受け入れて、その状態を上手く扱うために波動関数にラベルを付け加えたり、スピン変数（スピン座標）を付け加えたりしている、ということにして計算方法をインストールして頂ければと思います。^[なお、この「角運動量」というのは何か抽象的なものではなく、マクロな系の「回転」と同じものであることは確かなようです。参考：たとえばhttps://www.px.tsukuba.ac.jp/~onoda/ssh/node12.html]

わけわからんですね。しかし、すでに「物質は実は波」みたいな最大限にわけわからん仮定を受け入れて量子力学を使っているわけで、このような考え方も受け入れながら徐々に慣れていくしかないのかと思っています。また具体的にスピン関数を扱う中で、わけわからんなりに道具としては使えるようになります。
:::

新たにスピンを引数として加えた波動関数$\varphi(\boldsymbol{r},\sigma)$は、以下のように設定することで上手く扱うことができるようになります。なお、以降その関数を「**スピン軌道関数**」等と呼び、座標$\boldsymbol{r}$の関数を「軌道部分」や「軌道関数」と呼ぶことと対応し、スピンの関数を「スピン部分」「スピン関数」などと呼びます

まず「スピン変数」$\sigma$は、上向き、下向きに対応した2つの値を持ち、それを

$$
\sigma = \pm \frac{1}{2}
$$

と書いたり、

$$
\sigma = \uparrow, \downarrow
$$

と書いたりします。^[多分上向き・下向きの状態に紐づいていれば数値が何であっても問題ないと思うのですが（$\uparrow,\downarrow$とかありますし）、数値で表す際は$\pm\frac{1}{2}$とするのが慣例のようです。記事の後半でいずれ触れるつもりですが、球面調和関数$Y_l^m$の$z$方向角運動量演算子$\hat{L}_z$
の固有値、つまり固有状態$\varphi_{nlm}=R_{nl}Y_l^m$（$n,l$は任意）の、$z$方向の角運動量）が、量子数$m$を用いて$\hbar m$となることに対応して、$z$方向のスピン角運動量演算子$\hat{s}_z$の固有値が$\pm\frac{\hbar}{2}$である場合は、量子数$m_s=\pm\frac{1}{2}$に対応します。この$\pm\frac{1}{2}$を「スピン変数」としている感じだと思います。]

スピン関数の変数として扱う場合は、数値$\pm 1/2$を使い、状態のラベルとして扱う場合は$\uparrow,\downarrow$を使うことが多いです。本稿でも原則そのように使い分けるつもりです。


上向き状態に対応する「スピン関数」を$\alpha(\sigma)$, 下向き状態に対応するスピン関数を$\beta(\sigma)$と置き、それぞれスピン座標$\sigma=\pm 1/2$に対する値を

$$
\alpha\left(\frac{1}{2}\right) = 1, \alpha\left(-\frac{1}{2}\right)  = 0,\\
\beta\left(\frac{1}{2}\right) = 0, \beta\left(-\frac{1}{2}\right)  = 1,\\
$$

と置くと、

軌道部分が$\varphi(\boldsymbol{r})$でスピン状態が上向きの状態を

$$
\varphi_\uparrow(\boldsymbol{r},\sigma) = \varphi(\boldsymbol{r})\alpha(\sigma)
$$

軌道部分が$\varphi(\boldsymbol{r})$でスピン状態が下向きの状態を


$$
\varphi_\downarrow(\boldsymbol{r},\sigma)= \varphi(\boldsymbol{r})\beta(\sigma)
$$

と書けます。実際、$\varphi_\uparrow(\boldsymbol{r},\sigma)$、$\varphi_\downarrow(\boldsymbol{r},\sigma)$はそれぞれ、

$$
\begin{align*}
\varphi_\uparrow\left(\boldsymbol{r},\frac{1}{2}\right) &= \varphi(\boldsymbol{r})\alpha\left(\frac{1}{2}\right) = \varphi(\boldsymbol{r}),\\
\varphi_\uparrow\left(\boldsymbol{r},-\frac{1}{2}\right) &= \varphi(\boldsymbol{r})\alpha\left(-\frac{1}{2}\right) = 0,\\
\varphi_\downarrow\left(\boldsymbol{r},\frac{1}{2}\right) &= \varphi(\boldsymbol{r})\beta\left(\frac{1}{2}\right) = 0,\\
\varphi_\uparrow\left(\boldsymbol{r},-\frac{1}{2}\right) &= \varphi(\boldsymbol{r})\alpha\left(-\frac{1}{2}\right) = \varphi(\boldsymbol{r}),\\
\end{align*}
$$

なので、それぞれ$\varphi_\uparrow(\boldsymbol{r},\sigma)$は軌道部分が$\varphi(\boldsymbol{r})$で上向きスピンをもつ状態（略記した表現だと$\varphi_{\uparrow}$）、$\varphi_\downarrow(\boldsymbol{r},\sigma)$は軌道部分が$\varphi(\boldsymbol{r})$で下向きスピンをもつ状態（略記した表現だと$\varphi_{\downarrow}$）を上手く表現できていることがわかります。

### 内積と規格直交性

スピン関数同士の内積は$\alpha(\sigma)$を例にすると

$$
\begin{align*}
\int\alpha^*(\sigma)\alpha(\sigma)d\sigma &=\sum_\sigma\alpha^*(\sigma)\alpha(\sigma)\\
 &=\alpha^*\left(\frac{1}{2}\right)\alpha\left(\frac{1}{2}\right) + \alpha^*\left(-\frac{1}{2}\right)\alpha\left(-\frac{1}{2}\right)
\end{align*}
$$

のように定義され、$\alpha(\sigma), \beta(\sigma)$は以下の規格直交性を満たすものとします。

$$
\int\alpha^*(\sigma)\alpha(\sigma)d\sigma = 
\int\beta^*(\sigma)\beta(\sigma)d\sigma = 1,\\
\int\alpha^*(\sigma)\beta(\sigma)d\sigma = \int\beta^*(\sigma)\alpha(\sigma) d\sigma= 0
$$

またスピン関数の一般の状態として重ね合わせ状態^[注意点としては、この状態を観測した場合、スピンは中途半端な値を持つわけではなく、あくまで観測値は「上向き」または「下向き」のどちらかとなります。「$d_+$、$d_-$」の係数が意味するのは「同じ状態をたくさん集めて観測した場合に上向き（下向き）が観測される割合」です。
とはいえこのあたりはあまり深入りせず、「重ね合わせ状態もあるよ」という程度に認識しておいていただければ。]を考えることもでき、$|d_+|^2 + |d_-|^2 = 1$の係数を用いて

$$
\gamma(\sigma) = d_+\alpha(\sigma) + d_-\beta(\sigma)
$$

と書き、この関数も規格化されています：

$$
\int\gamma^*(\sigma)\gamma(\sigma)d\sigma = |d_+|^2 + |d_-|^2 = 1
$$

## ここまでのまとめ

というわけで今までの電子の状態についての考え方：「シュレーディンガー方程式を解いて得られた固有状態（やその重ね合わせ状態）によって記述される」に対して、以下の設定が新たに追加されることとなりました。

- 電子の状態は、位置座標に関する固有状態に加えて、スピンに関する2つの固有状態を持つ

- その状態をラベルを用いて$\varphi_{n\sigma}$と書いたり、スピン関数を用いて$\varphi(\boldsymbol{r})\gamma(\sigma)$と書いたりする

これにより前章で整理した多電子系の波動関数が満たすべき「反対称性」の表現が以下のように修正されることになります。

# スピンを考慮した多電子系の波動関数

## 前章の再掲

前章で、多電子系の波動関数について、以下のように整理しました。

::: message
波動関数は全電子の座標変数の関数

$$
\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots)
$$

によって表される。
:::

::: message
多電子系の波動関数は「反対称性」すなわち粒子の入れ替え

$$
\Phi(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N)
\rightarrow
\Phi(\boldsymbol{r}_2, \boldsymbol{r}_1,\cdots,\boldsymbol{r}_N)
$$

に対して波動関数に係数$-1$がかかる
:::


その上で、特に多粒子系を扱う場合でも粒子間の相互作用を扱わず、ハミルトニアンが一粒子の演算子の和

$$
\mathcal{H} \simeq \sum_i\hat{H}_i
$$

で書けた場合の波動関数の表現について以下のように整理しました。

::: message
一体のシュレディンガー方程式

$$
\hat{H}_i\phi_i(\boldsymbol{r}_i)  = \epsilon_i\phi_i(\boldsymbol{r}_i)
$$

を解いて固有関数を求め、得られた固有関数の積をスレーター行列式を用いて

$$
\begin{align*}
\Phi_\mathrm{F}(\boldsymbol{r}_1, \boldsymbol{r}_2,\cdots,\boldsymbol{r}_N) &= 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_a(\boldsymbol{r}_1) & \varphi_b(\boldsymbol{r}_1) & \cdots & \varphi_n(\boldsymbol{r}_1)\\
\varphi_a(\boldsymbol{r}_2) & \varphi_b(\boldsymbol{r}_2) & \cdots & \varphi_n(\boldsymbol{r}_2)\\
& \cdots & \cdots\\
\varphi_a(\boldsymbol{r}_N) & \varphi_b(\boldsymbol{r}_N) & \cdots & \varphi_n(\boldsymbol{r}_N)
\end{vmatrix}
\end{align*}
$$

と反対称化することで多粒子系の波動関数が得られる
:::

また重要な性質として、「パウリの排他原理」
::: message
スレーター行列式に含まれる1粒子固有関数$\varphi_a, \varphi_b,\cdots,\varphi_c,\cdots$はすべて異なる関数でなければならない
:::
を得ました。

この結果、多電子系の波動関数の基底状態は、

::: message
一粒子のシュレーディンガー方程式を解いて得られた固有エネルギーが小さい順から固有関数を選び、スレーター行列式に"詰めて”行くことで得られる
:::

と結論付けられました。

## スピンを考慮した多電子系の波動関数

本章で導入したスピンを考慮することで、上記の多電子系の波動関数が以下のように書き換えられます。

### 多電子系の波動関数
::: message
波動関数は全電子の座標変数$\boldsymbol{r}_i$とスピン変数$\sigma_i$を合わせた変数$\tau_i = (\boldsymbol{r}_i,\sigma_i)$の関数

$$
\Phi(\tau_1, \tau_2,\cdots)
$$

によって表される。
:::
### 反対称性
::: message
多電子系の波動関数は「反対称性」すなわち粒子の**位置座標とスピン座標**の入れ替え

$$
\Phi(\tau_1, \tau_2,\cdots)
\rightarrow
\Phi(\tau_2, \tau_1,\cdots)
$$

に対して波動関数に係数$-1$がかかる
:::

### 一体ハミルトニアンの和からなるシュレーディンガー方程式
ハミルトニアンが一粒子の演算子の和、**特にスピン変数に作用する演算子を含まない**、

$$
\mathcal{H} \simeq \sum_i\hat{H}_i,\\
\hat{H}_i = \hat{H}(\boldsymbol{r}_i,\nabla_i)
$$

で書けた場合（第二式は一粒子演算子が作用する粒子に対して形を変えない、の意）、

::: message
一体の（スピンを考慮しない）シュレディンガー方程式

$$
\hat{H}(\boldsymbol{r},\nabla)\phi(\boldsymbol{r})  = \epsilon\phi(\boldsymbol{r})
$$

を解いて**軌道部分の**固有関数$\varphi_a, \varphi_b\cdots$を求め、 **それらの上向き/下向きスピン状態$\varphi_{1\uparrow}, \varphi_{1\downarrow},\varphi_{2\uparrow}, \varphi_{2\downarrow}\cdots$、または別の表現をするとスピン関数$\alpha(\sigma)$、$\beta(\sigma)$の積** を一粒子固有関数とする。

**得られたスピン軌道関数を、軌道関数の量子数（例えば水素様原子中の電子の場合だと$(n,l,m)$）とスピン状態のラベル$\uparrow,\downarrow$をまとめて、ギリシャ文字$\lambda,\mu,\xi$などと表し^[本稿では可能な限り、軌道関数の固有状態のラベルをアルファベット$a,b,c\cdots$で書き、スピン状態も含めた固有状態のラベルをギリシャ文字$\lambda,\mu,\xi\cdots$で書くことにする。とはいえ明確な場合以外は都度どちらを意味しているか明記するようにする。]、位置座標$\boldsymbol{r}$とスピン座標$\sigma$を合わせて$\tau$と表す。**

スピン軌道関数$\varphi_\lambda(\tau_1), \varphi_\mu(\tau_2),\cdots,\varphi_\xi(\tau_N)$をスレーター行列式を用いて

$$
\begin{align*}
\Phi_\mathrm{F}(\tau_1, \tau_2,\cdots,\tau_N) &= 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_\lambda(\tau_1) & \varphi_\mu(\tau_1) & \cdots & \varphi_\xi(\tau_1)\\
\varphi_\lambda(\tau_2) & \varphi_\mu(\tau_2) & \cdots & \varphi_\xi(\tau_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\tau_N) & \varphi_\mu(\tau_N) & \cdots & \varphi_\xi(\tau_N)
\end{vmatrix}
\end{align*}
$$

と反対称化することで多粒子系の波動関数が得られる

スレーター行列式は、スピン関数との積で書く表現や、スピン関数$\alpha$との積の場合（上向きスピン状態の場合）は軌道関数のまま$\varphi_n$、スピン関数$\beta$との積（下向きスピン状態の場合）は軌道関数の上にバーをつけた$\overline{\varphi}_n$で書く表現などがある。また、含まれる一粒子固有状態の種類が指定できれば良いので、

$$
\begin{align*}
\Phi_\mathrm{F}(\tau_1, \tau_2,\cdots,\tau_N) &= 
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_\lambda(\tau_1) & \varphi_\mu(\tau_1) & \cdots & \varphi_\xi(\tau_1)\\
\varphi_\lambda(\tau_2) & \varphi_\mu(\tau_2) & \cdots & \varphi_\xi(\tau_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\tau_N) & \varphi_\mu(\tau_N) & \cdots & \varphi_\xi(\tau_N)
\end{vmatrix}\\
&\equiv \left| \varphi_a \overline{\varphi}_b\cdots \overline{\varphi}_n\right|
\end{align*}
$$

のように書くことも多い。本稿でもスペースとタイピング節約のために以降この略記法を多用する。（この略記法は複数の教科書で使われているが、一般的に通用するかどうかは保証されていない）

なお、これも後に出てくるので先に注意しておくと、**スレーター行列式の線形結合もまた、**

$$
\begin{align*}
\Phi_\mathrm{F}^{12}(\tau_1, \tau_2,\cdots,\tau_N)&\equiv C_1\Phi_\mathrm{F}^1(\tau_1, \tau_2,\cdots,\tau_N) + 
C_2 \Phi_\mathrm{F}^2(\tau_1, \tau_2,\cdots,\tau_N)\\
\Rightarrow \Phi_\mathrm{F}^{12}(\tau_2, \tau_1,\cdots,\tau_N) &=  C_1\Phi_\mathrm{F}^1(\tau_2, \tau_1,\cdots,\tau_N) + 
C_2 \Phi_\mathrm{F}^2(\tau_2, \tau_1,\cdots,\tau_N)\\
&=  - C_1\Phi_\mathrm{F}^1(\tau_1, \tau_2,\cdots,\tau_N) + 
-C_2 \Phi_\mathrm{F}^2(\tau_1, \tau_2,\cdots,\tau_N)\\
&= - \Phi_\mathrm{F}^{12}(\tau_1, \tau_2,\cdots,\tau_N)
\end{align*}
$$

**のように反対称性を満たし、多電子系の波動関数の解となり得る。**

:::

### パウリの排他原理

::: message
スレーター行列式に含まれる1粒子固有関数は、**軌道部分の固有関数
$\varphi_n(\boldsymbol{r})$ 1つに対して、上向きスピン状態$\varphi_{n\uparrow}$または$\varphi_n(\boldsymbol{r}\alpha(\sigma)$、下向きスピン状態$\varphi_{n\downarrow}$または$\varphi_n(\boldsymbol{r}\beta(\sigma)$の2つのみ（2つまで）** を取ることができる。
:::

なお、スピンも考慮に入れることで、パウリの排他原理の座標バージョンのような結果「同じスピン状態を持つ粒子は同じ座標を取らない」という効果も生まれます。（証明は追記予定）


### 多電子系の波動関数の基底状態

::: message
一粒子のシュレーディンガー方程式を解いて得られた固有エネルギーが小さい順から固有関数を選び、**スピン上向き状態、下向き状態をひとつづつ、計2つを選んで**スレーター行列式に"詰めて”行くことで多電子系の波動関数の基底状態が得られる。
※ここではハミルトニアン演算子にスピンの項が含まれない、すなわちスピンの向きによってエネルギーが変わらない場合を想定している。例えば磁場中の電子や、スピン軌道相互作用などハミルトニアン中にスピンに作用する演算子を含む場合はこの限りではない。そしておそらく本稿ではそういう状態は扱わない。
:::

これでようやく、スピンまで考慮に入れた場合の多電子系の波動関数の条件、具体的な問題の考え方が揃いました。

## 例題：N個の自由電子の波動関数

最後に改めて前章で扱った例題を、スピンまで考慮して解きなおしてみます。

N個の自由電子のハミルトニアン

$$
\mathcal{H} = \sum_{i=1}^N\left( \frac{-\hbar^2}{2m}\nabla_i^2 \right)
$$

を一体のハミルトニアンに分解して一体のシュレディンガー方程式

$$
\frac{-\hbar^2}{2m}\nabla^2 \phi(\boldsymbol{r}) = \epsilon\phi(\boldsymbol{r})
$$

を解き、固有関数と固有エネルギーを以下のように得るところまでは同じです。

$$
\varphi_{\boldsymbol{k}} = \frac{1}{\sqrt{V}}e^{i\boldsymbol{k}\cdot\boldsymbol{r}},\\
\epsilon_{\boldsymbol{k}} = C,\\
\boldsymbol{k} = (k_x, k_y, k_z), k_i = \frac{2\pi}{L}n_i,\\
(i = x,y,z, n_i = 0,\pm1,\pm2,\cdots)
$$

ここで、最低エネルギー状態（基底状態）の波動関数は、**固有エネルギーが小さい順から上向きスピン状態、下向きスピン状態をそれぞれひとつづつ**選んだスレーター行列式

$$
\Phi_\mathrm{F}(\tau_1, \tau_2,\cdots,\tau_N) = 
\left|
\varphi_{\boldsymbol{k}_0} \overline{\varphi}_{\boldsymbol{k}_0} \cdots  \overline{\varphi}_{\boldsymbol{k}_{N/2}}
\right|
$$


で書かれます。

固有エネルギーは一粒子シュレディンガー方程式を解いて得られた固有エネルギーを小さい順に足し合わせた

$$
E = \epsilon_{\boldsymbol{k}_0} + \epsilon_{\boldsymbol{k}_1} + \cdots + \epsilon_{\boldsymbol{k}_{N/2}}
$$

となります。（エネルギーが小さい順から軌道関数を並べて、一つの軌道状態当たりスピン上向き、下向きの2通りが取れるので、$N$粒子を考えた場合は$N/2$番目の固有エネルギーまでの和となります。）

### フェルミ波数・フェルミエネルギー

（追記予定）


# おわりに

ここまでで、なるべく角運動量などに触れずにスピンを導入し、スピンを考慮した多電子系の波動関数の性質を整理しました。おそらく「Tight-bindingモデルの理解」という目標だけを目指すのであれば、これで最低限の内容はかけたかと思います。

ここで整理した内容をもとに、次章ではいよいよ多電子原子中の電子の波動関数について扱っていきます。