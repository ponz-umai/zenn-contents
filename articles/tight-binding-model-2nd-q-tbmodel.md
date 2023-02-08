---
title: "Tight-binding modelの第二量子化表示"
emoji: "🎉"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: false
---

# はじめに

ついに最終章として、前章でまとめた第二量子化表示を用いて、本稿のゴールであるTight-bindingモデルの第二量子化表示を導出していきます。
最初はサクッと終わらせるつもりが、書き始めると気になることが出てきたり、やっぱりわからないことが出てきたり、そもそも記事を書くこと自体に想定以上に時間がかかったりと、だいぶ長くなってしまいましたが、ようやくこの章で一旦完結です。

第二量子化表示のTight-bindingモデルは、おそらく物性に関連する研究をしていたら見ない日は無いだろうと思われる、

$$
\mathcal{H} = \sum_{ i,j ,\sigma } -t_{ij}a^\dagger_{i\sigma }a_{j\sigma }  
$$

のように「飛び移り積分」（や"Hopping Integral", "Transfer Integral"等）と呼ばれる謎の定数$t$と、サイト$i,j$に電子を作ったり消したりする謎の演算子「生成・消滅演算子」によって構成される謎の式です。

おおよそ上式を眺めながら、

- サイト$i,j$に電子を「生成する」「消滅させる」って何？そんなことできるの？
- あとそもそもサイト$i,j$も何？どういう理屈で空間が離散化されてんの？（実はされてません）
- 場の演算子ってやつ？でもそっちは連続変数の位置座標$\boldsymbol{r}$に電子を生成するって書いてあったけど？それを離散化したやつなの？（違います）
- それでその係数の飛び移り積分ってのも何？
- ていうか全部何？

等と浮かんでは消える疑問から目を背けつつ、目の前に積みあがる研究タスクを消化するために「反交換関係」とやらを駆使しながら過ごす日々をお過ごしの方が、きっとこのような怪しげな個人ブログに行きついてきているものと思います。（そうですよね？それともそんな疑問に頭を悩ませていたのは私だけ？）

そんな日々に終止符を打つため、本章ではこれまでまとめてきた内容をもとに、上記のような第二量子化表示のTight-bindingモデルについて整理していきます。

# 前章の振り返り

初めに今回使う道具として、[前章](https://zenn.dev/ponzumai/articles/tight-binding-model-2nd-q)でまとめた第二量子化表示の基本的な考え方をこちらにまず再掲します。

まず、多体電子系のシュレーディンガー方程式を

$$
\mathcal{H}\Phi(\tau_1,\tau_2,\cdots,\tau_N)
=
E
\Phi(\tau_1,\tau_2,\cdots,\tau_N)
$$

において、任意の多体の波動関数$\Phi(\tau_1,\tau_2,\cdots,\tau_N)$を完全正規直交系$\{\varphi_\nu\}$で構成されるスレーター行列式の線形結合として

$$
\Phi(\tau_1,\tau_2,\cdots,\tau_N)
=
\sum_{\lambda,\mu,\cdots,\xi}C_{\lambda,\mu,\cdots\xi}
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

と展開すると考えると、「完全正規直交系$\{\varphi_\nu\}$にハミルトニアンを作用させるとどうなるか？」を考え、何かしらの演算子を対応させることでハミルトニアンを書き換えることができます。

特に多体のハミルトニアン$\mathcal{H}$が1電子のハミルトニアン$\hat{H}$の和で

$$
\mathcal{H} = \sum_i\hat{H}_i
$$

と書かれている場合、ハミルトニアンの作用は

$$
H_{\delta\kappa}\equiv \int \varphi_{\delta}^*(\tau)\hat{H}\varphi_{\kappa}(\tau)d\tau

$$

として、



$$
\begin{align*}
&\mathcal{H}\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|\\

&=\sum_i \hat{H}_i
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|\\
&=
\left|
\left(\hat{H}\varphi_\lambda\right)\varphi_\mu\cdots\varphi_\xi
    \right|

+
\left|
\varphi_\lambda\left(\hat{H}\varphi_\mu\right)\cdots\varphi_\xi
    \right|

+\cdots

+
\left|
\varphi_\lambda\varphi_\mu\cdots\left(\hat{H}\varphi_\xi\right)
    \right|\\

&=
\sum_{\delta}H_{\delta\lambda}\left|
    \varphi_{\delta}\varphi_\mu\cdots\varphi_\xi
    \right|
    +
    \sum_{\delta}H_{\delta\mu}\left|
    \varphi_{\lambda}\varphi_\delta\cdots\varphi_\xi
    \right|
    +\cdots
     +
    \sum_{\delta}H_{\delta\xi}\left|
    \varphi_{\lambda}\varphi_\mu\cdots\varphi_\delta
    \right|

\end{align*}
$$

となります。（2体相互作用の場合は省略）

そこで、以下のように生成演算子$\hat{a}_\delta^\dagger$と消滅演算子$\hat{a}_\kappa$を定義することで：

:::message

### 生成演算子の定義

スレーター行列式の一番手前に状態$\varphi_\delta$を付け加える：

$$
\hat{a}_\delta^\dagger
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|

=
\left|
\varphi_\delta\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

またはスレーター行列式に対応するラベルの1電子状態があった場合生成演算子$\hat{a}_\delta^\dagger$を作用させた結果は$0$となる：

$$
\hat{a}_\delta^\dagger
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\delta\cdots \varphi_\xi
    \right|

=
\left|
\varphi_\delta\varphi_\lambda\varphi_\mu\cdots\varphi_\delta\cdots \varphi_\xi
    \right|
=0.
$$

### 消滅演算子の定義


スレーター行列式の先頭にある1電子状態$\varphi_\kappa$を消す、先頭になければ先頭まで交換して移動してから消す：

$$
\hat{a}_\kappa|\varphi_\kappa\cdots \cdots|=|\cdots \cdots|
$$

または、

$$
\hat{a}_\kappa|(n個の状態)\varphi_\kappa \cdots|=
\hat{a}_\kappa(-1)^n|\varphi_\kappa(n個の状態) \cdots|=(-1)^n|\cdots \cdots|
$$

（この定義は$\hat{a}_\kappa\hat{a}_\kappa^\dagger = 1$とも見れる（多分））

また、スレーター行列式に消滅演算子$\hat{a}_\kappa$の示す状態$\varphi_\kappa$が含まれない場合は、消滅演算子を作用させた結果は$0$を返す：

$$
\hat{a}_\kappa|\cdots(\varphi_\kappa 無し) \cdots|=0.
$$

:::


多体のハミルトニアンのスレーター行列式への作用は二つの演算子を用いて

$$
\mathcal{H}\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|=\sum_i \hat{H}_i
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
    
=    
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

と書くことができます。

また、上記のように定義した生成・消滅演算子は、以下のような交換の規則「反交換関係」を満たします：

:::message
## 生成消滅演算子の満たす反交換関係

$$
\hat{a}_\kappa\hat{a}_\zeta + \hat{a}_\zeta\hat{a}_\kappa  
= 
\hat{a}_\delta^\dagger\hat{a}_\eta^\dagger + \hat{a}_\eta^\dagger\hat{a}_\delta^\dagger=0\\
\hat{a}_\delta^\dagger\hat{a}_\kappa + \hat{a}_\kappa\hat{a}_\delta^\dagger=
\delta_{\delta\kappa}
$$

:::


最後に、関数の変換

$$
\varphi_\nu = \sum_\eta c_{\eta\nu}\phi_\eta = \sum_\eta \braket{ \phi_\eta|\varphi_\nu}\phi_\eta
$$

逆変換

$$
\phi_\eta = \sum_\nu \tilde{c}_{\nu\eta}\phi_\eta= \sum_\nu \braket{ \varphi_\nu|\phi_\eta}\phi_\eta
$$


で結びつけられた2つの関数系$\{\varphi_\nu\}$、$\{\phi_\eta\}$を考え、それぞれ$\varphi_\nu$に対応した演算子を$\hat{a}_\nu,\hat{a}^\dagger_\nu$、変更後の$\phi_\eta$に対応した演算子を$\hat{b}_\eta, \hat{b}_\eta^\dagger$と置くと、両者は


:::message 
## 生成・消滅演算子の変換
$$
\hat{a}^\dagger_\nu
=
\sum_\eta c_{\eta\nu}b_\eta^\dagger
=
\sum_\eta \braket{ \phi_\eta|\varphi_\nu}\hat{b}_\eta^\dagger.
$$

$$
\hat{a}_\nu = \sum_\eta c_{\eta,\nu}^* \hat{b}_\eta = \sum_\eta \braket{\phi_\eta|\varphi_\nu} \hat{b}_\eta
$$

:::

として変換できます。

ここで、特に1体のハミルトニアンがスピンに関する項を含まない場合、第二量子化表示のハミルトニアンは


$$
\begin{align*}
\mathcal{H} &= \sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\\

&\equiv
\sum_{d,k, \gamma_d=\alpha,\beta, \gamma_k=\alpha,\beta}
\delta_{\gamma_d,\gamma_k}
\int \varphi_{d}^*(\boldsymbol{r})\hat{H}\varphi_{k}(\boldsymbol{r})d\boldsymbol{r}
\hat{a}^\dagger_{d,\gamma_d} \hat{a}_{k,\gamma_k}\\

&\equiv
\sum_{d,k, \gamma=\alpha,\beta}
H_{dk}
\hat{a}^\dagger_{d,\gamma} \hat{a}_{k,\gamma}

\end{align*}
$$

と、スピン部分については同じスピン状態を取り、軌道関数でハミルトニアンを挟んだ積分を係数とする形で書くことができます。

ここで2つめの式で出てきた$\hat{a}^\dagger_{d,\gamma}$は、スピン軌道関数のうち軌道関数部分が$\varphi_d(\boldsymbol{r})$、スピン関数部分がそれぞれ$\gamma=\alpha$なら$\alpha(\sigma)$となるスピン軌道関数$\varphi_\delta(\tau)=\varphi_d(\boldsymbol{r})\alpha(\sigma)$、$\gamma=\beta$なら$\beta(\sigma)$となるスピン軌道関数$\varphi_\delta(\tau)=\varphi_d(\boldsymbol{r})\beta(\sigma)$を生成する演算子です。（消滅演算子$\hat{a}_{k,\gamma}$も同様。）

なお、スピン部分のラベルは、$\sigma$を用いて$\sum_\sigma\hat{a}_{d,\sigma}^\dagger\hat{a}_{k,\sigma}$と書くことも多いです。

以上が前章で得られた第二量子化表示の基本的な内容の振り返りです。本章では、この知識を使いながらTight-bindingモデルの第二量子化表示を導出していきます。

# Tight-bindingモデルの第二量子化表示

以上のように整理してみると、第二量子化表示のために考えるべきことは

- (1) 多体電子系のハミルトニアンを設定する
- (2) スレーター行列式を展開するための何かしらの完全正規直交関数系を決める
- (3) 完全正規直交関数系でハミルトニアンを挟んで積分したハミルトニアン行列を求める

ことだとわかります。この3つが決まれば、後はこれまでに求めた手続きに従うことで第二量子化表示ができそうです。なお、2つ目の「完全正規直交関数系」を「基底」、「基底関数」などと呼び、特定の関数を選ぶことを「基底を決める」などと呼びます。この辺はベクトルっぽい用語が示すように「関数のベクトル表現」とかいう話と密接に関わる話なのですが、本稿では~~ギリ触れなくても説明できそうな気がしたので~~触れていません。が、「基底」という用語は字数も少なく済むので適宜使用します。^[皆さんご存じの内容とは思いますが、「関数のベクトル表現」「Hilbert空間」「Fock空間」などのキーワードで調べたり、適当な量子力学の教科書を読んだりするとちゃんと勉強できると思います。]

また、3番目のハミルトニアンの積分を求める際に、第一量子化で考えたようなTight-binding近似を取り入れていくことになります。

早速一つ目から始めていきましょう。

## (1) 多体電子系のハミルトニアンを設定する：固体中の多体電子系のハミルトニアン


さて、まずはハミルトニアンの設定から始めます。

これも基本的にこれまでの章、特に第一量子化のTight-bindingモデル[（前編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-1)[（後編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-2)で扱ってきた内容の復習です。

まずは上記の章で述べたように、大きな目標としては固体、すなわち周期的な格子ポテンシャルが存在するような空間の中に、多数の電子がいるような物理系における電子の振る舞い、すなわち波動関数を求めたいわけです。
しかしながら、ありのまますべて正確に書こうとすると、固体内の相互作用はこんな感じで

![](/images/tb/many-ele-atm-2.png)

多数の原子核の相互作用を受けながら、多数の電子がそれぞれに相互作用をしあっている状態です。わけわからんですね。ハミルトニアンもそれに対応して、

（原子核（や分子）からの相互作用ポテンシャル）＋（原子核同士の相互作用ポテンシャル）＋（電子間の相互作用ポテンシャル）＋（その他スピンに関する相互作用やらなんやかんや）

のように到底解けない形になってしまうわけです。

### ポテンシャルの近似

そこで、まず第一の近似として、

- 原子核（格子点）は静止しているものと近似する（Born–Oppenheimer近似）
- 同じ原子核内の多電子の相互作用は（例えば[多電子原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)の章で考えたような手法で）1体ポテンシャルとして近似し、結晶の周期と同じ周期で局所的な孤立原子のポテンシャルとして扱う
- 価電子間のCoulomb相互作用（やその他の相互作用）はひとまず無視する

のようにして、**お互いに相互作用しない電子が、周期的に並んだ（静止した）局所ポテンシャル（孤立原子ポテンシャル）を受けながら運動する**モデルを考えることにするわけです。これを1体近似などと呼びます。

これにより、孤立原子ポテンシャルを$V(\boldsymbol{r})$、結晶の周期を表す基本格子ベクトルを$\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3$と書き、各格子点の位置ベクトルを$\boldsymbol{R} = n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3, n_i = 0,\pm 1, \pm 2\cdots$と書いて、多体のハミルトニアンを

$$
\mathcal{H} \simeq \sum_i\left\{

\frac{-\hbar^2}{2m}\nabla_i{}^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}_i - \boldsymbol{R})

    \right\}
    \equiv \sum_i\hat{H}^{\rm c}_i
$$

と近似します。これでハミルトニアンを前章で整理した1体演算子の総和の形に書くことができました。

## (2) スレーター行列式を展開する完全正規直交関数系を決める

ハミルトニアンが1体演算子の形で書けたので、次はスレーター行列式を展開するための完全正規直交関数系を決めましょう。どんな関数系がこの条件を満たすでしょうか？

結論から言うと格子点$\boldsymbol{R}$を中心とする原子軌道関数（孤立原子ポテンシャルの固有関数）$\{\phi_m(\boldsymbol{r} - \boldsymbol{R})\}$で展開していくことになります。

### 周期的境界条件と平面波展開

関数系を考える前に、境界条件を設定します。ここではよく使われる（というか「表面状態を見たい」等特別な理由がない限りほぼ使われている）Born–von Karmanの周期的境界条件と呼ばれる境界条件を設定します。これは周期的な境界条件ではありますが、その周期を十分大きく取ることで結局とても（無限に）大きい結晶を考えるような条件です。具体的には、先ほど周期的なポテンシャルにおいて設定した格子の周期を表す基本格子ベクトルと「十分大きな整数」$N_i$を用いて、波動関数が

$$
\varphi(\boldsymbol{r} + N_1\boldsymbol{a}_1) = \varphi(\boldsymbol{r} + N_2\boldsymbol{a}_2) = \varphi(\boldsymbol{r} + N_3\boldsymbol{a}_3) = \varphi(\boldsymbol{r} ) 
$$

を満たすと設定します。

このように設定することで、波動関数はいわば「$N$倍の基本格子ベクトル」とでもいうべき周期を持つ周期関数と考えられて、この時[格子ベクトルと逆格子ベクトル](https://zenn.dev/ponzumai/articles/tight-binding-model-lattice-vec)の章で見たように、以下のような波数ベクトル$\boldsymbol{q}$

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
\varphi(\boldsymbol{r}) = \frac{1}{\sqrt{V}}\sum_{\boldsymbol{q}}c_{\boldsymbol{q}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}
$$

ここで$\frac{1}{\sqrt{V}}$は規格化のための係数で、$V$は結晶の単位胞の体積$v_c = \boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)$の$N=N_1N_2N_3$倍の定数$V = v_cN_1N_2N_3$です。

つまり、「完全正規直交関数系」第一候補は平面波$\{\frac{1}{\sqrt{V}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}\}$です。^[Fourier展開を用いて任意の（正確には区分的なめらかな）周期関数を表すことができるという完全性は証明せずに使います]

とはいえ愚直にこの平面波で展開しても別にいいことはないわけで、もう少し考えていきます。

### 逆格子ベクトル

なお、この後使うのでついでに「逆格子ベクトル」について定義しておきます。詳細は[格子ベクトルと逆格子ベクトル](https://zenn.dev/ponzumai/articles/tight-binding-model-lattice-vec)をご確認ください。

先ほど1体ポテンシャルを設定した際に用いた結晶の周期を表す基本格子ベクトルを$\boldsymbol{a}_i$として、その任意の整数倍の和「格子ベクトル」を

$$
\begin{align*}
\boldsymbol{R} &= n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3\\
\end{align*}
$$

としたとき、格子ベクトルの周期を持つ関数$V(\boldsymbol{r} + \boldsymbol{R}) = V(\boldsymbol{r})$のFourier級数展開は以下のように書けます。


$$
V(\boldsymbol{r}) = \sum_{\boldsymbol{K}}C_{\boldsymbol{K}}e^{i\boldsymbol{K}\cdot\boldsymbol{r}}\\
C_{\boldsymbol{K}}= \frac{1}{v_c}\int_{v_c} V(\boldsymbol{r})e^{-i\boldsymbol{K}\cdot\boldsymbol{r}}d\boldsymbol{r}
$$

この時の波数$\boldsymbol{K}$を逆格子ベクトルと呼び、$\boldsymbol{K}$は以下のように定義されます：

$$
\begin{align*}
    \boldsymbol{K} 

&=
 m_1\frac{2\pi\boldsymbol{a}_2\times\boldsymbol{a}_3}{\boldsymbol{a}_1\cdot(\boldsymbol{a}_2\times\boldsymbol{a}_3)}
+

 m_2\frac{2\pi\boldsymbol{a}_3\times\boldsymbol{a}_1}{\boldsymbol{a}_2\cdot(\boldsymbol{a}_3\times\boldsymbol{a}_1)}
+
 m_3\frac{2\pi\boldsymbol{a}_1\times\boldsymbol{a}_2}{\boldsymbol{a}_3\cdot(\boldsymbol{a}_1\times\boldsymbol{a}_2)}\\

 &\equiv
 m_1\boldsymbol{b}_1 + m_2\boldsymbol{b}_2 + m_3\boldsymbol{b}_3
\end{align*}
$$


### Blochの定理とBloch関数とWannier関数

次に、[固体（結晶）中の電子状態](https://zenn.dev/ponzumai/articles/tight-binding-model-electrons-in-solids)で整理したBlochの定理を思い出していきます。

#### Blochの定理（の一部）

今考えているような結晶の周期を持つ（1体）ハミルトニアン

$$
\hat{H}^{\rm c}
 =

\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})
$$

の固有関数は、上記の平面波$\{\frac{1}{\sqrt{V}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}\}$のセットから、逆格子ベクトルの最小周期、特に（第一）ブリルアンゾーンと呼ばれる領域の波数$\boldsymbol{k}$を選び、その波数$\boldsymbol{k}$と先ほど定義した逆格子ベクトル$\boldsymbol{K}$だけ異なる波数だけを抜き出した線形結合

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \frac{1}{\sqrt{V}}\sum_{\boldsymbol{K}}c_{n,\boldsymbol{k} - \boldsymbol{K}}
e^{i(\boldsymbol{k} - \boldsymbol{K} )\cdot \boldsymbol{r}}
$$

によって書くことができます。このように固有関数は、組み合わせる波数を代表する波数$\boldsymbol{k}$（この波数を「Bloch波数」と呼びます）と、（無限個の）逆格子ベクトルだけ異なる平面波の組み合わせから得られる（無限個の）固有値に対応したラベル$n$によって指定され、Bloch関数などと呼ばれます。
Bloch関数に対応する固有値を、ラベル$n,\boldsymbol{k}$を用いて$\varepsilon_{n,\boldsymbol{k}}$と書き、固有値、固有関数（Bloch関数）は$\hat{H}^{\rm c}\varphi_{n,\boldsymbol{k}} = \varepsilon_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}$を満たします。固有値のラベル$n$は、固有値（固有エネルギー）が小さい順に番号を振ることにしておきます。

#### 完全正規直交基底その2：Bloch関数

さて、前置きが長くなりましたが、上記のようにして構成されるBloch関数は正規直交性$\int\varphi^*_{n',\boldsymbol{k}'}(\boldsymbol{r})\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r} = \delta_{n',n}\delta_{\boldsymbol{k}',\boldsymbol{k}}$を満たしまた、完全正規直交性を満たす平面波のセット$\{\frac{1}{\sqrt{V}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}\}$の線形結合、いわば組み換えなわけなので任意の関数を表すことができる完全性も満たします。

:::details 一応証明しておきます

（本当はこんなことしなくても、エルミート演算子の固有関数系は完全性を満たすという定理があるようなのですが、なんか大変そうなのでまあそれは一回見なかったことにして。。。）（あと、以下のように考えておくとこの後のWannier関数の完全性も示せますし）

Bloch関数の定義

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \frac{1}{\sqrt{V}}\sum_{\boldsymbol{K}}c_{n,\boldsymbol{k} - \boldsymbol{K}}
e^{i(\boldsymbol{k} - \boldsymbol{K} )\cdot \boldsymbol{r}}
$$

と、平面波の規格直交性から、

$$
c_{n,\boldsymbol{k} - \boldsymbol{K}} = \int_V e^{-i(\boldsymbol{k}-\boldsymbol{K})\cdot \boldsymbol{r}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r}
$$

です。一方平面波からBloch関数への逆変換を

$$
e^{i(\boldsymbol{k} - \boldsymbol{K} )\cdot \boldsymbol{r}} = \frac{1}{\sqrt{V}}\sum_n d_{n,\boldsymbol{k}-\boldsymbol{K}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})
$$

と書くと、先ほどと同様にBloch関数の規格直交性から係数は

$$
d_{n,\boldsymbol{k}-\boldsymbol{K}} = \int_V \varphi_{n,\boldsymbol{k}}^*e^{i(\boldsymbol{k} - \boldsymbol{K} )\cdot \boldsymbol{r}} 
d\boldsymbol{r}

= c_{n,\boldsymbol{k}-\boldsymbol{K}}^*
$$

とわかります。

さて、Fourier級数の完全性から、任意の波動関数は

$$
\psi(\boldsymbol{r}) = \sum_{\boldsymbol{q}}A_{\boldsymbol{q}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}
$$

とあらわせますが、ここで波数の和をBloch波数と逆格子ベクトルで分けて、先ほどの平面波からBloch関数への変換を代入して、

$$
\begin{align*}
\psi(\boldsymbol{r}) &= \sum_{\boldsymbol{q}}A_{\boldsymbol{q}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}\\

&= 
\sum_{\boldsymbol{k}}\sum_{\boldsymbol{K}}A_{\boldsymbol{k}-\boldsymbol{K}}e^{i(\boldsymbol{k}-\boldsymbol{K})\cdot\boldsymbol{r}}\\

&=
\sum_{\boldsymbol{k}}\sum_{\boldsymbol{K}}A_{\boldsymbol{k}-\boldsymbol{K}}\frac{1}{\sqrt{V}}
\sum_n c_{n,\boldsymbol{k}-\boldsymbol{K}}^*\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})
\end{align*}
$$

なので、$\sum_{\boldsymbol{K}}\frac{1}{\sqrt{V}}A_{\boldsymbol{k}-\boldsymbol{K}}c_{n,\boldsymbol{k}-\boldsymbol{K}}^*\equiv B_{n,\boldsymbol{k}}$等と置きなおせば任意の関数をBloch関数を用いて

$$
\psi(\boldsymbol{r}) = \sum_{n,\boldsymbol{k}}B_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})
$$

と展開できることが示せました。（この証明合ってるのかちょっと自信ないですが。。。なんか間違ってることに気づいたら直します）

:::

というわけで、第二の完全正規直交関数系候補として、Bloch関数$\{\varphi_{n,\boldsymbol{k}}\}$が挙げられます。

なお、Bloch関数は1体ハミルトニアンの固有状態になっているので、これを基底に選べばハミルトニアンの第二量子化表示は冒頭で書いたハミルトニアンにスピンの項がない場合を用いて

$$
\begin{align*}
\mathcal{H} 
&=
 \sum_{\sigma=\uparrow,\downarrow}\sum_{n',\boldsymbol{k}'}\sum_{n,\boldsymbol{k}}\int\varphi_{n',\boldsymbol{k}'}^*\hat{H}^{\rm c}\varphi_{n,\boldsymbol{k}}d\boldsymbol{r}\hat{a}^\dagger_{n',\boldsymbol{k}',\sigma}\hat{a}_{n,\boldsymbol{k},\sigma}\\

&=
\sum_{\sigma=\uparrow,\downarrow}\sum_{n',\boldsymbol{k}'}\sum_{n,\boldsymbol{k}}
\delta_{n'n}\delta_{\boldsymbol{k},\boldsymbol{k}'}\varepsilon_{n,\boldsymbol{k}}\hat{a}^\dagger_{n',\boldsymbol{k}',\sigma}\hat{a}_{n,\boldsymbol{k},\sigma}
\\

&=
\sum_{\sigma=\uparrow,\downarrow}\sum_{n,\boldsymbol{k}}
\varepsilon_{n,\boldsymbol{k}}\hat{a}^\dagger_{n,\boldsymbol{k},\sigma}\hat{a}_{n,\boldsymbol{k},\sigma}

\end{align*}
$$

と、自明な固有状態を持つ形になります。これを「対角化された」形などと呼ぶことが多いです。Tight-bindingモデルに行き着く前に対角化表示が出てきてしまいましたが、結局Bloch関数は「平面波を何かしらの係数$c_{n,\boldsymbol{k}-\boldsymbol{K}}$で線形結合したもの」という情報しかなく、このままでは固有値も求められないわけで、固有値を求めるためにも具体的に計算可能な関数を用意したいわけです。

#### 完全正規直交関数：Wannier関数

ここで、[Bloch関数の局在関数ーWannier関数ーを用いた展開](https://zenn.dev/ponzumai/articles/tight-binding-model-wannier-func)の章で見たように、Bloch関数はBloch波数$\boldsymbol{k}$の周期性

$$
\varphi_{n,\boldsymbol{k} + \boldsymbol{K}}(\boldsymbol{r}) = \varphi_{n, \boldsymbol{k}}(\boldsymbol{r})
$$

あるいは、$N_i\rightarrow\infty$とすれば$\boldsymbol{k}$は連続変数となることを踏まえてBloch関数を$\varphi_n(\boldsymbol{r}, \boldsymbol{k})$と書いて

$$
\varphi_n(\boldsymbol{r},\boldsymbol{k} + \boldsymbol{K}) = \varphi_n(\boldsymbol{r}, \boldsymbol{k})
$$

を満たします。

このことから、波数ベクトル$\boldsymbol{k}$の関数$\varphi_n(\boldsymbol{r}, \boldsymbol{k})$は、「逆格子ベクトルの逆格子ベクトル」を波数として持つ**逆格子空間の**平面波を用いてFourier級数展開できますが、[逆格子ベクトルの章](https://zenn.dev/ponzumai/articles/tight-binding-model-reciprocal-lattice)で見た通り「逆格子ベクトルの逆格子ベクトル」は、実空間の格子ベクトル、つまり

$$
\begin{align*}
\boldsymbol{R} &= n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3\\
\end{align*}
$$

と一致するので、結局逆格子空間の平面波

$$
e^{i\boldsymbol{R}\cdot\boldsymbol{k}}
$$


を用いてFourier展開（平面波展開）

$$
\varphi_n(\boldsymbol{r},\boldsymbol{k}) = \sum_{\boldsymbol{R}}C_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{R}\cdot\boldsymbol{k}},
\\

C_{n,\boldsymbol{R}}(\boldsymbol{r})= \frac{1}{v_{BZ}}\int_{BZ} \varphi_n(\boldsymbol{r},\boldsymbol{k})e^{-i\boldsymbol{R}\cdot\boldsymbol{k}}d\boldsymbol{k}
$$


できるのでした。ここで2式目の積分は逆格子空間の単位胞、すなわちブリルアンゾーンを積分範囲として取ります。また$v_{BZ}$はブリルアンゾーンの体積$\boldsymbol{b}_1\cdot(\boldsymbol{b}_2\times\boldsymbol{b}_3)$です。

そして、上記のように展開したBloch関数において、波数$\boldsymbol{k}$についての展開を考えていたことから、展開係数$C_{n,\boldsymbol{R}}$もまた、$\boldsymbol{r}$の関数となっており電子の観測確立に対応した波動関数的なものであると考えられます。このように定義された$C_{n,\boldsymbol{R}}(\boldsymbol{r})$を$w_{n,\boldsymbol{R}}(\boldsymbol{r})$と書き、この関数を提唱者の名前を取ってWannier関数と呼びます。改めて書くと、

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \varphi_n(\boldsymbol{r},\boldsymbol{k}) = \sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}},\\
w_{n,\boldsymbol{R}}(\boldsymbol{r})
=
\frac{1}{v_{BZ}}\int_{BZ} \varphi_n(\boldsymbol{r},\boldsymbol{k})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}d\boldsymbol{k}
$$


です。Wannier関数は実空間上の位置座標$\boldsymbol{r}$の関数であり、かつ実空間上の格子点の座標でラベルされるややこしい見た目をしています。


またここで「Fourier逆変換」部分$w_{n,\boldsymbol{R}}(\boldsymbol{r})=\frac{1}{v_{BZ}}\int_{BZ} \varphi_n(\boldsymbol{r},\boldsymbol{k})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}d\boldsymbol{k}$はWannier関数の定義とみることもできますが、$\boldsymbol{k}$を離散変数と考えて積分を和で書き換えることもあります。

具体的には、この時実空間の単位胞の体積を$v_{cell}$、周期的境界条件の周期を$N = N_1N_2N_3$とし、その1周期分（つまり固体全体）の体積を$V = Nv_{cell}$として

$$
\int_{BZ} d\boldsymbol{k}\rightarrow\frac{(2\pi)^3}{V}\sum_{\boldsymbol{k}\in BZ} = \frac{(2\pi)^3}{Nv_{cell}}\sum_{\boldsymbol{k}\in BZ} 
$$

となります。ここで$\sum_{\boldsymbol{k}\in BZ}$は（第1）ブリルアンゾーン内の$\boldsymbol{k}$について和を取ることを示します。さらにブリルアンゾーンの体積$v_{BZ}$が実空間の単位胞の体積$v_{cell}$を用いて

$$
v_{BZ} = \frac{2\pi}{v_{cell}}
$$

となることを用いて、Wannier関数の定義を

$$
w_{n,\boldsymbol{R}}(\boldsymbol{r})
=
\frac{1}{N}\sum_{\boldsymbol{k}\in BZ} \varphi_{n\boldsymbol{k}}(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
$$

と書くこともできます。（この式を利用する方が多いです）

::: message alert
文献によっては$1/N$を2式に振り分けて、

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \frac{1}{\sqrt{N}}\sum_{\boldsymbol{R}}w_{n,\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}},\\
w_{n,\boldsymbol{R}}(\boldsymbol{r})
=
\frac{1}{\sqrt{N}}\sum_{\boldsymbol{k}} \varphi_n(\boldsymbol{r},\boldsymbol{k})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
$$

と定義することも多いです。

:::

このようにして定義されたWannier関数は、以下の性質を満たします。
- 並進性：$w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_n(\boldsymbol{r}-\boldsymbol{R})$
- 正規直交性：$\int_V w_{n'}^*(\boldsymbol{r}-\boldsymbol{R}')w_n(\boldsymbol{r} - \boldsymbol{R})dr =\delta_{n,n'}\delta_{\boldsymbol{R},\boldsymbol{R}'}$
- 局在性：ラベル＝格子点の座標$\boldsymbol{R}$を中心として局在した関数


先ほど「Wannier関数は実空間上の位置座標$\boldsymbol{r}$の関数であり、かつ実空間上の格子点の座標でラベルされるややこしい見た目をしています」と書きましたが、ラベル$\boldsymbol{R}$は実空間上の関数の「中心」または「平行移動する分」を意味しているのですね。

そして、このようにして定義されたWannier関数

$$
\begin{align*}
w_{n,\boldsymbol{R}}(\boldsymbol{r})&=

\frac{1}{v_{BZ}}\int_{BZ} \varphi_n(\boldsymbol{r},\boldsymbol{k})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}d\boldsymbol{k}\\
&=
\frac{1}{N}\sum_{\boldsymbol{k}\in BZ} \varphi_{n\boldsymbol{k}}(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
\end{align*}
$$

もまた、完全正規直交関数系であるBloch関数$\{\varphi_{n\boldsymbol{k}}\}$の変換・逆変換ができるため、先ほど平面波$\rightarrow$Bloch関数で示したように完全性を満たす基底になり得ます。

とはいえそもそものBloch関数の具体的な関数系が指定されなければ、Wannier関数やWannier関数によるハミルトニアンの積分も計算できません。次へ進みましょう。


#### 完全正規直交関数：原子軌道関数と重なり積分の近似

さて、前置きが長くなりましたが、ようやく具体的な関数系を導入します。

Wannier関数の「局在性」に注目し、原点に局在するWannier関数を、同じく原点を中心とする完全正規直交関数系である「原子軌道関数」を用いて展開しよう、と考えるわけです。

ここで原子軌道関数を以下のように、孤立原子ポテンシャル1つだけがある系の固有関数として定義します。（孤立原子ポテンシャルは固体中のポテンシャルを$\hat{H}^{\rm c}=\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})$とを1体近似で近似した際、周期的に並んでいるとしたポテンシャル$V(\boldsymbol{r})$です。）



すなわち、孤立原子のハミルトニアンを$\hat{H}^{\rm a}=-(\hbar^2/2m)\nabla^2 + V(\boldsymbol{r})$として、

$$
\hat{H}^{\rm a} \phi_m(\boldsymbol{r}) = \varepsilon_m^{\rm a}\phi_m(\boldsymbol{r}). 
$$

ここで固有関数$\phi_m(\boldsymbol{r})$を「原子軌道」や「原子軌道関数」、固有値$\varepsilon_m^{\rm a}$を「原子準位」等と呼ぶことにします。

**今、孤立原子のポテンシャルは（例えば[多電子原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-many-electron-atom)の章で考えたようなHF近似や、さらに進んだ何らかの手法によって）具体的に分かっていると考えます。そしてその固有関数系も、頑張って具体的に得られたものと仮定します。すなわち、この先出てくる$\{\phi_m\}$は、今までのBloch関数やWannier関数とは異なり、具体的に（解析的にでも、数値的にでも）関数系が分かっている関数として考えていきます。**

このように定義した原子軌道関数を用いて、先ほどBloch関数を展開した「格子点$\boldsymbol{R}$に中心を持つ」Wannier関数$w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_{n}(\boldsymbol{r}-\boldsymbol{R})$を、「格子点$\boldsymbol{R}$に中心を持つ」孤立原子の固有関数$\phi_m(\boldsymbol{r}-\boldsymbol{R})$と展開係数$b_m$を用いて、

$$
w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_n(\boldsymbol{r}-\boldsymbol{R}) = \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R}) 
$$

と展開します。ここで、Wannier関数の並進性を考えると任意の$\boldsymbol{R}$ラベルの$w_{n,\boldsymbol{R}}$は、別のラベル$\boldsymbol{R}'$を平行移動しただけの関数ですので、展開係数$b_m^n$は$\boldsymbol{R}$によりません。ただ原子軌道関数が$\boldsymbol{R}$だけ平行移動されるだけとなります。
また展開係数は規格化条件$\sum_m|b_m^n|^2 = 1$を満たすものとします。

さらに略記法として原子準位と「中心の格子点の座標」で指定される原子軌道関数を

$$
\phi_m(\boldsymbol{r}-\boldsymbol{R})\equiv \phi_{m,\boldsymbol{R}}(\boldsymbol{r})
$$

と定義しておきます。

以上踏まえると、ここでも完全正規直交関数系Wannier関数から原子軌道関数への変換が定義されたので、Bloch関数やWannier関数の場合と同様に任意の関数を関数系$\{\phi_{m,\boldsymbol{R}}\}$で展開することができ、このセットをスレーター行列式を展開する基底に選ぶことができそうです。

ただし問題は、関数系$\{\phi_{m,\boldsymbol{R}}\}$は直交性を満たしていません。

どういうことかというと、同じ中心を持つ原子軌道同士はエルミート演算子である孤立原子のハミルトニアンの固有状態なので直交性$\int\phi_{m'}^*(\boldsymbol{r})\phi_m(\boldsymbol{r})d\boldsymbol{r}=\delta_{m'm}$を満たしますが、異なる中心の原子軌道同士が直交する保証はなく、

$$
\int\phi_{m',\boldsymbol{R}'}^*(\boldsymbol{r})\phi_{m,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}=

\int\phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}

\neq \delta_{m'm}\delta_{\boldsymbol{R}'\boldsymbol{R}}
$$

と直交性を満たしていません。とはいえ完全性は満たし、また規格化はされているので、イメージとしては斜交ベクトルみたいな感じなのかと思います（多分。これは適当なことを書いています。）

この場合、第二量子化のレシピに沿って考えて行く中で波動関数にハミルトニアンを作用させた結果を

$$
\hat{H}\phi_{m,\boldsymbol{R}} = \sum_{m'\boldsymbol{R}'}c_{m',\boldsymbol{R}'}^{m,\boldsymbol{R}}\phi_{m'\boldsymbol{R}'}
$$

と展開するところまでは同じようにできますが、その展開係数がハミルトニアン行列の形になりません。つまり

$$
c_{m'\boldsymbol{R}'}^{m,\boldsymbol{R}}\neq
\int \phi_{m'\boldsymbol{R}'}^*(\boldsymbol{r})\hat{H}\phi_{m,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}
$$

とできません。これは困ってしまいます。

世の中には中心が異なる原子軌道が直交性を満たすように原子軌道関数系を定義する方法もあるっぽいのですが^[または[格子ベクトルと逆格子ベクトル](https://zenn.dev/ponzumai/articles/tight-binding-model-lattice-vec)の場合の類推で考えると、何か「双対関数」とでも呼ぶべき関数$\tilde{\phi}_{m,\boldsymbol{R}}(\boldsymbol{r}): \int\tilde{\phi}_{m',\boldsymbol{R}'}^*(\boldsymbol{r})\phi_{m,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}=\delta_{m'm}\delta_{\boldsymbol{R}'\boldsymbol{R}}$を考えれば展開係数を$c_{m'\boldsymbol{R}'}^{m,\boldsymbol{R}}=\int \tilde{\phi}_{m'\boldsymbol{R}'}^*(\boldsymbol{r})\hat{H}\phi_{m,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}$と定義することができるかもしれないが、そんなことどこでも見たことないしどうやってそんな関数を作るのかもよくわからない。]、一旦そういうややこしいことは考えず、以下のように原子軌道間の**重なり積分**を定義し

$$
s_{\boldsymbol{R}',\boldsymbol{R}}^{m',m} = \int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\phi_m(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}
$$

異なるラベル間の重なり積分をゼロと近似すれば（このような近似は[重なり積分（Overlap Integral）の物理的意味（の妄想）](https://zenn.dev/ponzumai/articles/tight-binding-model-overlap-int)を考えることで、まあそれなりに妥当であろうと結論付けたのでした）：

$$
s_{\boldsymbol{R}',\boldsymbol{R}}^{m',m} = \int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\phi_m(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}
=\delta_{m'm}\delta_{\boldsymbol{R}'\boldsymbol{R}}
$$

原子準位$m$と関数の中心を表す格子点の座標$\boldsymbol{R}$で指定される原子軌道関数系は規格直交性を満たすようになります。


$$
\int\phi_{m',\boldsymbol{R}'}^*(\boldsymbol{r})\phi_{m,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}=

\int\phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}

\neq \delta_{m'm}\delta_{\boldsymbol{R}'\boldsymbol{R}}
$$

これでスレーター行列式を展開する基底（完全規格直交関数系）として、具体的な関数$\{\phi_{m,\boldsymbol{R}}\}$を採用することができました。

先述の通り原子軌道関数は具体的な関数を準備することができるので、これで結晶のハミルトニアン行列を計算することができ、第二量子化表示を具体的に書くことができるようになりました。

あともう少し考えることがあるのですが、ここで一回現状をまとめておきます。

第二量子化の考え方に沿って整理していくと、今各格子点を中心とする原子軌道関数$\{\phi_m(\boldsymbol{r} - \boldsymbol{R})\}$または$\boldsymbol{R}$もラベルとして書いて$\{\phi_{m,\boldsymbol{R}}(\boldsymbol{r})\}$を軌道関数部分の完全正規直交基底として選びました。

この時この関数系とスピン関数$\gamma=\alpha,\beta$どちらかの積を

$$
\phi_{m,\boldsymbol{R}, \uparrow}(\boldsymbol{r})\equiv
\phi_{m,\boldsymbol{R}}(\boldsymbol{r})\alpha(\sigma),\\

\phi_{m,\boldsymbol{R}, \downarrow}(\boldsymbol{r})\equiv
\phi_{m,\boldsymbol{R}}(\boldsymbol{r})\beta(\sigma)
$$

と書き、スピン状態を指定するラベル$\uparrow,\downarrow$をまとめて$\gamma$と書くことにして（$\phi_{m,\boldsymbol{R}, \gamma}$）、これらを並べて構成したスレーター行列式$|\phi_{m,\boldsymbol{R},\gamma}, \phi_{m',\boldsymbol{R}',\gamma'},\cdots\phi_{m'''',\boldsymbol{R}'''',\gamma''''}|$で多電子状態の波動関数を

$$
\Phi(\tau_1, \tau_2, \cdots, \tau_N) = \sum C_{m,m'\cdots;\boldsymbol{R},\boldsymbol{R}'\cdots;\gamma,\gamma'\cdots}
\left|
    \phi_{m,\boldsymbol{R},\gamma}, \phi_{m',\boldsymbol{R}',\gamma'},\cdots\phi_{m'''',\boldsymbol{R}'''',\gamma''''}
    \right|
$$

と展開すると考えることに対応しています。

この時第二量子化のレシピを利用して、特に1体のハミルトニアン$\frac{-\hbar^2}{2m}\nabla_i{}^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}_i - \boldsymbol{R})$にスピンを含む項がないため、


$$
\begin{align*}
\mathcal{H} &= \sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\\

&=
\sum_{m',\boldsymbol{R}';m,\boldsymbol{R} \gamma=\uparrow,\downarrow}
\int \phi_{m',\boldsymbol{R}'}^*(\boldsymbol{r})\hat{H}\phi_{m,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}
\hat{a}^\dagger_{m',\boldsymbol{R}'\gamma} \hat{a}_{m,\boldsymbol{R}\gamma}\\

&\equiv
\sum_{m,\boldsymbol{R};m,\boldsymbol{R} \gamma=\uparrow,\downarrow}
\braket{m',\boldsymbol{R}'|\hat{H}|m,\boldsymbol{R}}
\hat{a}^\dagger_{m',\boldsymbol{R}'\gamma} \hat{a}_{m,\boldsymbol{R}\gamma}\\

\end{align*}
$$

と書け、あとは（どこかでよく見た）$\int \phi_{m',\boldsymbol{R}'}^*(\boldsymbol{r})\hat{H}\phi_{m,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}$の積分を計算するだけとなりました。

これでかなりゴールに近づいた感があります。



## (3) ハミルトニアン行列を求める/LCAO近似/Tight-binding近似

最後に、冒頭で述べた第二量子化3つのステップの3番目、(3) ハミルトニアン行列の計算に取り掛かっていきましょう。ここで、これまでの章で考えてきた様々なTight-binding近似を取り入れていくことにします。

### LCAO近似

初めに「固体中の電子の固有状態は少数の原子軌道関数の重ね合わせで表せる」と考えるLCAO近似（Linear Combination of Atomic Orbitals、つまり「原子軌道の線形結合」の頭文字を取ってこのように呼ばれます）を考えます。

これは[第一量子化のTight-bindingモデル（前編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-1)で考えたように、まず固体のハミルトニアン$\hat{H}^{\rm c}=\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})$の固有関数で、固有値$\varepsilon_{n,\boldsymbol{k}}$に対応する$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$が原子軌道関数で

$$
\varphi_{n,\boldsymbol{k}}(\boldsymbol{r}) = \sum_m b_m^n \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot \boldsymbol{R}}\phi_m(\boldsymbol{r} - \boldsymbol{R})
$$

と展開します。（Bloch関数をWannier関数で展開して、展開したWannier関数を原子軌道関数で展開したらこうなります。また、上式は関数$\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot \boldsymbol{R}}\phi_m(\boldsymbol{r} - \boldsymbol{R})$の係数$b_m^n$の線形結合とみることもできます。この$\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot \boldsymbol{R}}\phi_m(\boldsymbol{r} - \boldsymbol{R})$の形を「Bloch和」と呼び、Blochの定理を満たします。ただし固体のハミルトニアンの固有関数ではありません。）



ここで、原子軌道関数$\phi_m$の孤立原子のハミルトニアンに対する固有値（原子準位）を$\varepsilon_m^{\rm a}$として、原子準位$\varepsilon_l^{\rm a}\sim \varepsilon_{n,\boldsymbol{k}}$を持つ原子軌道に掛かる係数$b^n_l$と、原子準位$\varepsilon_L^{\rm a}\gg \varepsilon_{n,\boldsymbol{k}}$を持つ原子軌道に掛かる係数$b^n_L$との関係を

$$
\begin{align*}
\frac{b^n_L}{b^n_l} &\simeq

(\varepsilon_{n,k} - \varepsilon_l^{\rm a})\frac{-S_L^n }{V_l^n}

\simeq 0

\end{align*}
$$

と近似できることに対応しています。なお、

$$
S_l^n\equiv\sum_m 
\left[
   \sum_{\boldsymbol{R}\neq\boldsymbol{0}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} 
   \int \phi_l^*(\boldsymbol{r})
   \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
   \right]b_m^n
$$


$$
V_l^n\equiv

\sum_m 
\left[
   \sum_{\boldsymbol{R}}
   e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \int \phi_l^*(\boldsymbol{r})
   \left(
   \sum_{\boldsymbol{R}'\neq \boldsymbol{0}}
      V(\boldsymbol{r}-\boldsymbol{R}')
      \right)  \phi_m(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
      \right]b^n_m 
$$

です。この近似は「$V_l^n$が$S_l^n$に比べてそれなりに大きい場合」に成り立つ近似であると考えられます。それぞれの積分の意味を考えてみると、これはつまり「局所ポテンシャルがそれなりに大きい場合」と見ることができ、結局「Tight-bindingな場合」の近似に対応しているとわかります。

この近似により、以下のように考えることができます。

原子軌道関数のうちエネルギー的に近い原子準位を持つ原子軌道のラベルの集合

$$
Q = \{m_{q_1},m_{q_2},\cdots m_{q_{|Q|}}\}
$$

を考えます。（ラベルの付け方をスッキリさせるために「集合」等と書いてますが要はグループにして名前を付けただけです）ここで$Q$に含まれる原子軌道の数を$|Q|$と書いています

$Q$に属する原子軌道の原子準位を$\varepsilon_{m_{q_1}}^{\rm a}, \varepsilon_{m_{q_2}}^{\rm a},\cdots \varepsilon_{m_{q_{|Q|}}}^{\rm a}$として、これら原子準位と近い固有値を$\varepsilon_{n_q,\boldsymbol{k}} \simeq \varepsilon_{{q_1}}^{\rm a}, \cdots$と書きます。

以上のように考えると、$|Q|$個の原子軌道関数から作られたBloch和の重ね合わせから、$|Q|$個の固体の固有関数が以下のように**近似的に**展開されます
（ラベルの付け方のセンスが絶望的で申し訳ありません。何か改善案を思いついたらもう少し見やすくします）：
$$
\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{m \in Q}b_m^{n_{p_i}}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R}),\\
i = 1,2,\cdots |Q|.
$$

ここで、この展開によって得られる固有関数のラベルの集合を

$$
P = \{n_{p_1}, n_{p_2}, \cdots \} 
$$

と置きます。

この時詳細は[第一量子化のTight-bindingモデル（前編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-1)を見ていただくとして、近似的な固有関数$\varphi_{n_i,\boldsymbol{k}}$は行列の一般化固有値問題の解として得られるため、固有関数は自動的に「エネルギーが近い」と選んだ原子軌道の数$|Q|$と同じ数だけ得られます。つまり$|P|=|Q|$を満たします。

また「どれくらい近ければ同じ集合に入れるのか？」と疑問に思うかもしれません。実際の研究でどうやっているのかは分かりませんが、あくまで近似なので「実験で得られたエネルギーバンドを再現できるまで」ということになるかとは思います。例えば教科書には、「遷移金属のバンド構造を計算するときには、一般には$d$と$s$の両方の準位」^[アシュクロフト・マーミン上(I)p.241]である6準位を一つのエネルギーが近いグループとして考える等と書いてあります。
また[水素原子中の電子](https://zenn.dev/ponzumai/articles/tight-binding-model-hydrogen-atom)で書いたように、原子準位はその角度部分のラベル$s,p,d\cdots$ごとに$1,3,5,\cdots$個に縮退しているので、そのうち一つ縮退したの軌道群を考えるのが最も素朴な近似で、少し先取りして書くとこのように一つの準位（例えば$s$）のみとりだして得られるバンドを「$s$バンド」等と呼び、この時に定義される生成消滅演算子を「$s$電子の生成消滅演算子」等と呼んだりします。この辺はまた後で再度言及します。


話を戻して、上記のように少数の原子軌道で固体の固有関数を展開できると仮定した場合、その逆変換を以下のように定義しておきます。（ここで先述のように異なる格子点に中心を持つ原子軌道の重なり積分をゼロと近似しています）

$$
\phi_{m_{q_j}}(\boldsymbol{r} - \boldsymbol{R}) = \frac{1}{N}\sum_{n\in P}\tilde{b}_{m_{q_j}}^{n}\sum_{k\in BZ}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\varphi_{n,\boldsymbol{k}},\\

j = 1,2,\cdots|Q|
$$


### ハミルトニアン行列を求める：Tight-binding近似の適用


これで、残すはハミルトニアン行列要素

$$
\begin{align*}
\braket{m,\boldsymbol{R}'|\hat{H}|m,\boldsymbol{R}} &= \int \phi_{m',\boldsymbol{R}'}^*(\boldsymbol{r})\hat{H}\phi_{m,\boldsymbol{R}}(\boldsymbol{r})d\boldsymbol{r}\\

&=
\int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\hat{H}\phi_{m}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
\end{align*}
$$

の計算となりました。

この積分はすでに[飛び移り積分（Hopping Integral）の物理的意味・Wannier関数の従う方程式](https://zenn.dev/ponzumai/articles/tight-binding-model-hopping-int)の章で考えています。

違いとしては今回はLCAO近似すなわち固有関数を少数の原子軌道で展開できるという近似を考えているため、あるグループ$Q$に属する、格子点$\boldsymbol{R}$を中心に持つ原子軌道関数にハミルトニアンを作用させた関数

$$
\begin{align*}
\hat{H}^{\rm c} \phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R}) 
=
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}-\boldsymbol{R}'')
\right) \phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R}) 
\end{align*}
$$


は、以下のようにして**同じ「エネルギー的に近いグループ」に属する原子軌道**の、原子準位と格子点の位置についての線形結合で展開できます：
$$
\begin{align*}
\hat{H}^{\rm c}\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})&=

\sum_{m'\in Q,\boldsymbol{R}'}c_{m',\boldsymbol{R}'}^{m_{q_j},\boldsymbol{R}}\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}').
\end{align*}
$$

::: 証明

まず上式に「逆変換」$\phi_{m_{q_j}}(\boldsymbol{r} - \boldsymbol{R}) = \frac{1}{N}\sum_{n\in P}\tilde{b}_{m_{q_j}}^{n}\sum_{k\in BZ}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\varphi_{n,\boldsymbol{k}}$を代入すると、$\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})$は$\hat{H}^{\rm c}$の固有関数なので、


$$
\begin{align*}
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}-\boldsymbol{R}'')
\right) \phi_{m_{q_j}}(\boldsymbol{r}) 

&=
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}-\boldsymbol{R}'')
\right) 
\frac{1}{N}\sum_{n\in P,\boldsymbol{k}\in BZ}\tilde{b}_{m_{q_j}}^n
e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})
\\

&=
\frac{1}{N}\sum_{n\in P,\boldsymbol{k}}\varepsilon_{n,\boldsymbol{k}}\tilde{b}_{m_{q_j}}^ne^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})

\end{align*}
$$

となりますが、さらに今度は$\varphi_{n\boldsymbol{k}}$を原子軌道関数（Bloch和）で展開して、

$$
\begin{align*}
（上式右辺）&=
\frac{1}{N}\sum_{n\in P,\boldsymbol{k}}\varepsilon_{n,\boldsymbol{k}}\tilde{b}_{m_{q_j}}^n\sum_{m'\in Q} b_{m'}^n\sum_{\boldsymbol{R}'}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}'-\boldsymbol{R})} \phi_{m'}(\boldsymbol{r} - \boldsymbol{R}')
\end{align*}
$$

さらにBloch波数の関数$\varepsilon_{n,\boldsymbol{k}}$の^[Fourier展開（Bloch関数に加え、固有関数もBloch波数について逆格子ベクトルを周期とする周期関数なのでした。というわけで、Bloch関数と同様に、直接格子の格子ベクトルを用いてFourier展開できます）]

$$
\varepsilon_{n,\boldsymbol{k}} = \sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{R}},\\
\varepsilon_{n,\boldsymbol{R}} = \frac{1}{N}\sum_{\boldsymbol{k}}e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}\varepsilon_{n,\boldsymbol{k}} 
$$

を代入して、

$$
\begin{align*}
（上式右辺）&=
\frac{1}{N}\sum_{n\in P,\boldsymbol{k}}\sum_{\boldsymbol{R}''}e^{i\boldsymbol{k}\cdot\boldsymbol{R''}}\varepsilon_{n,\boldsymbol{R}''}\tilde{b}_{m_{q_j}}^n\sum_{m'\in Q} b_{m'}^n\sum_{\boldsymbol{R}'}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}'-\boldsymbol{R})} \phi_{m'}(\boldsymbol{r} - \boldsymbol{R}')



\end{align*}
$$

となります。ここで$\boldsymbol{k}$の総和が計算できるので、$\sum_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot(\boldsymbol{R}' + \boldsymbol{R}''-\boldsymbol{R})} = N\delta_{\boldsymbol{R}-\boldsymbol{R}',\boldsymbol{R}''}$を用いて最終的に、

$$
\begin{align*}
\hat{H}^{\rm c}\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})&=
\sum_{m'\in Q}
\sum_{\boldsymbol{R}'}

\left(\sum_{n\in P}\varepsilon_{n,\boldsymbol{R}-\boldsymbol{R}'}\tilde{b}_{m_{q_j}}^n b_{m'}^n
\right)\phi_{m'}(\boldsymbol{r} - \boldsymbol{R}')


\\
&\equiv
\sum_{m'\in Q,\boldsymbol{R}'}c_{m',\boldsymbol{R}'}^{m_{q_j},\boldsymbol{R}}\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}').
\end{align*}
$$

:::


今、重なり積分をゼロと置いているので、展開係数は


$$
\begin{align*}

c_{m',\boldsymbol{R}'}^{m_{q_j},\boldsymbol{R}}

&=
\left\{
\begin{array}{ll}
\int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\hat{H}\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}& (m'\in Q )\\
0 & (m'\notin Q)
\end{array}
\right.

\\

&\equiv

\delta_{m'\in Q}'
\int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\hat{H}\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}

\end{align*}
$$


であると考えることができます。なお、ここで拡張されたクロネッカーデルタとして、

$$
\delta_{m'\in Q}' \equiv
\left\{
\begin{array}{ll}
1& (m'\in Q )\\
0 & (m'\notin Q)
\end{array}
\right.
$$

を定義しました。多分一般的には使われない表現だと思うので注意してください。（ちょっと探したのですが今回のように「変数がある集合に含まれるときだけ値を持つ」みたいな場合の書き方をどうすれば良いのかわからなかったのでとりあえず上記のように書きましたが、正しい（？）記法がわかれば修正します）

最後に$m'\in Q$を満たす場合の積分

$$
\int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\hat{H}\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}

=
\int\phi^*_{m'}(\boldsymbol{r} - \boldsymbol{R}')\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}-\boldsymbol{R}'')
\right)\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
$$



の具体的な計算をします。と言ってもこれも[飛び移り積分（Hopping Integral）の物理的意味・Wannier関数の従う方程式](https://zenn.dev/ponzumai/articles/tight-binding-model-hopping-int)の章ですでに考えており、その結果を引用すると、


まず、孤立原子ハミルトニアン部分の固有関数であることを利用して、右側の原子軌道にハミルトニアンの孤立原子部分を作用させて固有値を抜き出し、

$$
\begin{align*}
（上式右辺）&=

\varepsilon_{m_{q_j}}^{\rm a}\int\phi_{m'}^*(\boldsymbol{r} + \boldsymbol{R}')\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d \boldsymbol{r} 
+ 
\int\phi_{m'}^*(\boldsymbol{r} - \boldsymbol{R}')
\left(
  \sum_{\boldsymbol{R}'' \neq \boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}'')
\right)
\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d \boldsymbol{r} \\

&=
\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}\delta_{\boldsymbol{R}'\boldsymbol{R}}+ 
\int\phi_{m'}^*(\boldsymbol{r} - \boldsymbol{R}')
\left(
  \sum_{\boldsymbol{R}'' \neq \boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}'')
\right)
\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d \boldsymbol{r} 

\end{align*}
$$


ここで$\boldsymbol{R}' = \boldsymbol{R}$とその他で分けると、

$\boldsymbol{R}' = \boldsymbol{R}$の場合：

$$
\begin{align*}
c_{m',\boldsymbol{R}}^{m,\boldsymbol{R}} &\simeq

\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}

+
\int\phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R})
\left(
  \sum_{\boldsymbol{R}'' \neq \boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}'')
\right)
\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d \boldsymbol{r} \\

&=
\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}

+\Delta\varepsilon_{m'{m_{q_j}}}
\end{align*}
$$

と、原子準位$\varepsilon_{m_{q_j}}^{\rm a}$と結晶場積分$\Delta\varepsilon_{m'{m_{q_j}}}$で書けることがわかります。

ここで、多くの場合結晶場積分は[第一量子化のTight-bindingモデル（後編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-2)の章でちらっと述べたようにラベルによらない定数と近似されます。（一応今は残しておきます）


$\boldsymbol{R}' \neq \boldsymbol{R}$の場合：

また$\boldsymbol{R}' \neq \boldsymbol{R}$の場合は、Tight-bindingモデルの第一量子化で行ったように「3中心積分」をゼロと置くと、


$$
\begin{align*}
c_{m',\boldsymbol{R}'\neq \boldsymbol{R}}^{{m_{q_j}},\boldsymbol{R}} 
&=
\int\phi_{m'}^*(\boldsymbol{r} - \boldsymbol{R}')
\left(
  \sum_{\boldsymbol{R}'' \neq \boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R}'')
\right)
\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d \boldsymbol{r} \\



&\simeq

\int\phi_{m'}^*(\boldsymbol{r} - \boldsymbol{R}')
  V(\boldsymbol{r}-\boldsymbol{R}')
\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d \boldsymbol{r} 
\\

&\equiv

-t_{(m',\boldsymbol{R}') \leftarrow ({m_{q_j}},\boldsymbol{R})}
\\
&=
-t_{(m',\boldsymbol{R}'-\boldsymbol{R}) \leftarrow ({m_{q_j}},\boldsymbol{0})}

\end{align*}
$$

ここで、飛び移り積分を

$$
-t_{(m,\boldsymbol{R}') \leftarrow (n,\boldsymbol{R})}
\equiv

\int
 \phi_m^*(\boldsymbol{r}-\boldsymbol{R}')
 
   V(\boldsymbol{r} - \boldsymbol{R}')
 
 \phi_n(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}

$$

と定義しました。なお、飛び移り積分はポテンシャルの周期性と周期的境界条件から、格子点の座標ベクトルの差$\boldsymbol{R}'-\boldsymbol{R}$のみに依存するため、$t_{(m',\boldsymbol{R}') \leftarrow ({m_{q_j}},\boldsymbol{R})}=t_{(m',\boldsymbol{R}'-\boldsymbol{R}) \leftarrow ({m_{q_j}},\boldsymbol{0})}$が成り立ちます。（なので「前後の座標ベクトルの差の関数」みたいにかければいいのですが、そうするとなんかわかりにくくなるので辞めました）


またラベル$(m,\boldsymbol{R}') \leftarrow (n,\boldsymbol{R})$は、[飛び移り積分（Hopping Integral）の物理的意味・Wannier関数の従う方程式](https://zenn.dev/ponzumai/articles/tight-binding-model-hopping-int)の章で~~妄想~~考察したように、「状態$m$,格子点$\boldsymbol{R}$にいる原子軌道が、微小時間後に状態$n$で格子点$\boldsymbol{R}'$に「飛び移っている」確率（に比例する量）と解釈できることに対応しています。

最後に、全ての飛び移り積分を考えるのではなく、**隣接格子**間の飛び移り積分のみが値を持つと考える近似を行います。

この時飛び移り積分は二つの格子点の相対ベクトルのみに依存するので、相対ベクトルだけで指定することにします。

すなわち、ある格子点$\boldsymbol{R}$から、必要に応じた分だけの（多くは最近接か、次近接くらいまでのようです）隣接格子への相対ベクトルを$\boldsymbol{N}$と置くと、ある格子点で状態$n$から、隣接格子ベクトルだけ離れた格子点$\boldsymbol{R}+\boldsymbol{N}$へ、状態$m$となる飛び移り積分を、

$$
\begin{align*}
-t_{(m,\boldsymbol{R}+\boldsymbol{N}) \leftarrow (n,\boldsymbol{R})}
&=

\int
 \phi_m^*(\boldsymbol{r}-\boldsymbol{N})
 
   V(\boldsymbol{r} - \boldsymbol{N})
 
 \phi_n(\boldsymbol{r})d\boldsymbol{r}\\

 &=

 \int
 \phi_m^*(\boldsymbol{r}-(\boldsymbol{R}+\boldsymbol{N}))
 
   V(\boldsymbol{r} - (\boldsymbol{R}+ \boldsymbol{N}))
 
 \phi_n(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}

\end{align*}
$$

と書けます。ここで、「隣接ベクトル」$\boldsymbol{N}$は1つではないことに注意してください。（例えば1次元格子をイメージしてもらうと、最低でも「右隣」と「左隣」の二方向を考える必要があります）


この時$\boldsymbol{N}$として$\boldsymbol{N}_1, \boldsymbol{N}_2\cdots \boldsymbol{N}_i\cdots$を考えるとして、$\boldsymbol{N}_i$のラベルの集合を$I$とおくと
^[隣接格子ベクトルのラベルに$I$がよく使われるので。でも何の頭文字化はよくわかっていません。]、「飛び移り後」の格子点の座標$\boldsymbol{R}'$と、飛び移り前の格子点の座標$\boldsymbol{R}$の差（相対ベクトル）が$\boldsymbol{N}_i, i\in I$のどれかに一致するとき以外ゼロになる、として

$$
\begin{align*}
c_{m',\boldsymbol{R}'\neq \boldsymbol{R}}^{{m_{q_j}},\boldsymbol{R}} 
&\equiv

-t_{(m',\boldsymbol{R}') \leftarrow ({m_{q_j}},\boldsymbol{R})}
\\

&\simeq

-t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m_{q_j}},\boldsymbol{R})}
\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}

\end{align*}
$$

と定義しておくと良さそうです。

以上より、ハミルトニアン行列は、$\hat{H}^{\rm c}\phi$の展開係数$c_{m',\boldsymbol{R}'\neq \boldsymbol{R}}^{{m_{q_j}},\boldsymbol{R}}$を、$\boldsymbol{R}=\boldsymbol{R}$の時にゼロ、それ以外で$1$となるように

$$
(1-\delta_{\boldsymbol{R}\boldsymbol{R}'})
$$

をかけて、

$$
\begin{align*}
c_{m',\boldsymbol{R}'\neq \boldsymbol{R}}^{{m_{q_j}},\boldsymbol{R}} 

&\equiv

-t_{(m',\boldsymbol{R}') \leftarrow ({m_{q_j}},\boldsymbol{R})}
(1-\delta_{\boldsymbol{R}\boldsymbol{R}'})

\end{align*}
$$

あるいは隣接格子間の飛び移り積分だけ考えると、さらに近似して

$$
\begin{align*}
c_{m',\boldsymbol{R}'\neq \boldsymbol{R}}^{{m_{q_j}},\boldsymbol{R}} 

&\equiv

-t_{(m',\boldsymbol{R}') \leftarrow ({m_{q_j}},\boldsymbol{R})}
(1-\delta_{\boldsymbol{R}\boldsymbol{R}'})\\

&\simeq

-t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m_{q_j}},\boldsymbol{R})}
\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}

\end{align*}
$$



とすると$\boldsymbol{R}=\boldsymbol{R}', \boldsymbol{R}\neq \boldsymbol{R}'$の時をまとめて以下のように

$$
\begin{align*}
&\braket{m'\boldsymbol{R}'|\hat{H}^{\rm c}|{m_{q_j}}\boldsymbol{R}}\\
&=\int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\hat{H}\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
\\
&=
\int\phi^*_{m'}(\boldsymbol{r} - \boldsymbol{R}')\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}-\boldsymbol{R}'')
\right)\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
\\

&\simeq
\left\{
\left(
\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}

+\Delta\varepsilon_{m'{m_{q_j}}}
\right)\delta_{\boldsymbol{R}'\boldsymbol{R}}


-t_{(m',\boldsymbol{R}') \leftarrow ({m_{q_j}},\boldsymbol{R})}(1-\delta_{\boldsymbol{R}'\boldsymbol{R}})
\right\}\delta_{m'\in Q}'\\

&\simeq
\left\{
\left(
\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}

+\Delta\varepsilon_{m'{m_{q_j}}}
\right)\delta_{\boldsymbol{R}'\boldsymbol{R}}

-t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m_{q_j}},\boldsymbol{R})}
\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}
\right\}\delta_{m'\in Q}'

\end{align*}
$$

と表せます。

## Tight-bindingモデルの第二量子化表示

以上で、第二量子化表示に必要な3ステップが以下のようにすべて準備できました。

### (1) 多体電子系のハミルトニアンを設定する：
:::message
孤立原子のポテンシャルが、固体内の格子点に周期的に並んでいると考え、かつ電子間のCoulomb相互作用を無視して

$$
\mathcal{H}

\simeq

\sum_i
  \left(
\frac{-\hbar^2}{2m}\nabla_i{}^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}_i-\boldsymbol{R}'')
\right) 

\equiv
\sum_i\hat{H}_i^{\rm c}
$$

と一体近似をしたハミルトニアンを考える。

:::
### (2) スレーター行列式を展開するための何かしらの完全正規直交関数系を決める：

:::message
固体のハミルトニアンを$\hat{H}^{\rm c}=\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})$と1体近似で近似した際、周期的に並んでいるとしたポテンシャル$V(\boldsymbol{r})$が「孤立して」ある場合のハミルトニアン（孤立原子のハミルトニアン）$\hat{H}^{\rm a}=-(\hbar^2/2m)\nabla^2 + V(\boldsymbol{r})$の固有関数：

$$
\hat{H}^{\rm a} \phi_m(\boldsymbol{r}) = \varepsilon_m^{\rm a}\phi_m(\boldsymbol{r}). 
$$

およびその格子ベクトル分平行移動した関数系

$$
\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

なお、LCAO近似・Tight-binding近似より、この原子軌道関数のうちエネルギー的に近い原子準位を持つ原子軌道のラベルの集合

$$
Q = \{m_{q_1},m_{q_2},\cdots m_{q_{|Q|}}\}
$$

を考え、$Q$に属する原子軌道の原子準位を$\varepsilon_{m_{q_1}}^{\rm a}, \varepsilon_{m_{q_2}}^{\rm a},\cdots \varepsilon_{m_{q_{|Q|}}}^{\rm a}$とし、これら原子準位と近い固有値をもつ固体のハミルトニアンの固有状態$\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r})$：$\hat{H}^{\rm c}\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n_{p_i},\boldsymbol{k}}\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r})$、$\varepsilon_{n_q,\boldsymbol{k}} \simeq \varepsilon_{{q_1}}^{\rm a}, \cdots$とする。

すると$|Q|$個の原子軌道関数から作られたBloch和の重ね合わせから、$|Q|$個の固体の固有関数を以下のように**近似的に**展開できる：

$$
\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{m \in Q}b_m^{n_{p_i}}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R}),\\
i = 1,2,\cdots |Q|.
$$

また、この展開によって得られる固有関数のラベルの集合を

$$
P = \{n_{p_1}, n_{p_2}, \cdots \} 
$$

と置き、$|P|=|Q|$を満たす。


:::
  
### (3) 完全正規直交関数系でハミルトニアンを挟んで積分したハミルトニアン行列を求める

:::message
上記の完全正規直交関数系と、固有関数の近似的な展開の関係を用いて、

$$
\begin{align*}
&\braket{m'\boldsymbol{R}'|\hat{H}^{\rm c}|{m_{q_j}}\boldsymbol{R}}\\
&=\int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\hat{H}\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
\\
&=
\int\phi^*_{m'}(\boldsymbol{r} - \boldsymbol{R}')\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}-\boldsymbol{R}'')
\right)\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
\\

&\simeq

\left\{
\left(
\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}

+\Delta\varepsilon_{m'{m_{q_j}}}
\right)\delta_{\boldsymbol{R}'\boldsymbol{R}}


-t_{(m',\boldsymbol{R}') \leftarrow ({m_{q_j}},\boldsymbol{R})}(1-\delta_{\boldsymbol{R}'\boldsymbol{R}})
\right\}\delta_{m'\in Q}'\\

&\simeq
\left\{
\left(
\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}

+\Delta\varepsilon_{m'{m_{q_j}}}
\right)\delta_{\boldsymbol{R}'\boldsymbol{R}}


-t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m_{q_j}},\boldsymbol{R})}
\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}
\right\}\delta_{m'\in Q}'
\end{align*}
$$

となる。
ここで、飛び移り積分を

$$
-t_{(m,\boldsymbol{R}') \leftarrow (n,\boldsymbol{R})}
\equiv

\int
 \phi_m^*(\boldsymbol{r}-\boldsymbol{R}')
 
   V(\boldsymbol{r} - \boldsymbol{R}')
 
 \phi_n(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}

$$

定義した。この時隣接格子への飛び移り積分を、隣接格子への相対ベクトルを$\boldsymbol{N}$と置いて

$$
\begin{align*}
-t_{(m,\boldsymbol{R}+\boldsymbol{N}) \leftarrow (n,\boldsymbol{R})}
&\equiv

\int
 \phi_m^*(\boldsymbol{r}-\boldsymbol{N})
 
   V(\boldsymbol{r} - \boldsymbol{N})
 
 \phi_n(\boldsymbol{r})d\boldsymbol{r}\\

 &=

 \int
 \phi_m^*(\boldsymbol{r}-(\boldsymbol{R}+\boldsymbol{N}))
 
   V(\boldsymbol{r} - (\boldsymbol{R}+ \boldsymbol{N}))
 
 \phi_n(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}

\end{align*}
$$

と書け、どちらも前後の格子点の座標の差のみに依存する。
:::

最後に、これらを用いて、生成消滅演算子を用いた第二量子化表示を書いていきましょう。

完全正規直交関数系に対応して、関数$\phi_{m,\boldsymbol{R}}(\boldsymbol{r})=\phi_{m}(\boldsymbol{r}-\boldsymbol{R})$とスピン関数$\gamma(\sigma)=\alpha(\sigma), \beta(\sigma)$の積

$$
\phi_{m,\boldsymbol{R}}(\boldsymbol{r})\gamma(\sigma)
\equiv
\phi_{m,\boldsymbol{R},\gamma}(\boldsymbol{r})
$$

をスレーター行列式に「付け加える」（既にあれば$0$を返す）生成演算子を$\hat{a}_{m,\boldsymbol{R},\gamma}^\dagger$：

$$
\hat{a}_{m,\boldsymbol{R},\gamma}^\dagger|\cdots|=|\phi_{m\boldsymbol{R},\gamma}\cdots|
$$

スレーター行列式の先頭から$\phi_{m,\boldsymbol{R},\gamma}(\boldsymbol{r})$を消す（無ければ$0$を返す、先頭ではない場所にあれば先頭まで移動してから消す）消滅演算子を$\hat{a}_{m,\boldsymbol{R},\gamma}$：

$$

\hat{a}_{m,\boldsymbol{R},\gamma}|\phi_{m,\boldsymbol{R},\gamma}\cdots|=|\cdots|

$$

をを定義すると、固体のハミルトニアンの完全正規直交関数系$\{\phi_{m,\boldsymbol{R}}\}$で展開された任意のスレーター行列式への作用、つまり第二量子化表示は、あるラベル$m$が属する集合を$Q(\ni m)$と書くことにして、さらに全ての$m$の集合を$Q_1, Q_2\cdots Q_i \cdots$と分割することにして、（こんな書き方が正しいかどうかはわからないのですが）

$$
\begin{align*}
\mathcal{H}

&\simeq

\sum_i
  \left(
\frac{-\hbar^2}{2m}\nabla_i{}^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}_i-\boldsymbol{R}'')
\right) 

\\

&=\sum_{\gamma=\uparrow,\downarrow}\sum_{m'\boldsymbol{R}'m\boldsymbol{R}}\braket{m'\boldsymbol{R}'|\hat{H}^{\rm c}|{m}\boldsymbol{R}}
\hat{a}_{m',\boldsymbol{R}',\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}
\\



&\simeq
\sum_{\gamma=\uparrow,\downarrow}\sum_{m'\boldsymbol{R}'m\boldsymbol{R}}
\left[
\left\{
\left(
\varepsilon_{m}^{\rm a}\delta_{{m}m'}

+\Delta\varepsilon_{m'{m}}
\right)\delta_{\boldsymbol{R}'\boldsymbol{R}}


-t_{(m',\boldsymbol{R}') \leftarrow ({m},\boldsymbol{R})}(1-\delta_{\boldsymbol{R}'\boldsymbol{R}})
\right\}\delta_{m'\in Q(\ni m)}'
\right]
\hat{a}_{m',\boldsymbol{R}',\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}\\

&=

\sum_{Q_i}\left\{
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    \sum_{m'\in Q_i}
    \left(
        \sum_{\boldsymbol{R}}
        
        \left(
    \varepsilon_{m}^{\rm a}\delta_{{m}m'}

    +\Delta\varepsilon_{m'{m}}
    \right)
    \hat{a}_{m',\boldsymbol{R},\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}


    \sum_{\boldsymbol{R}\neq\boldsymbol{R}'}
    -t_{(m',\boldsymbol{R}') \leftarrow ({m},\boldsymbol{R})}\hat{a}_{m',\boldsymbol{R}',\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}
    \right)
\right\}
\\

\end{align*}
$$

と、**「エネルギーが近いグループ」ごとに分かれた形**になります。^[これが特定のバンドの生成消滅演算子のみで第二量子化表示ができる理由です。私はこの部分が分からず、永らく頭を悩ませまくっていました。]

また先述のように、結晶場積分をラベル$m,m'$によらない定数と近似し、基準となるエネルギーを結晶場積分の値の分だけずらして式から消去して、かつ飛び移り積分を隣接格子だけ取ると近似すると、

$$
\begin{align*}
\mathcal{H}

&=

\sum_{Q_i}\left\{
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    
        \sum_{\boldsymbol{R}}
        
    \varepsilon_{m}^{\rm a}
    \hat{a}_{m,\boldsymbol{R},\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}

    +
    \sum_{\gamma=\uparrow,\downarrow}
\sum_{m\in Q_i}
\sum_{m'\in Q_i}

    \sum_{\boldsymbol{R}\neq\boldsymbol{R}'}
  -t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m},\boldsymbol{R})}
\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}
    
    \hat{a}_{m',\boldsymbol{R}',\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}
\right\}
\\

\end{align*}
$$

となります。ここで、$\sum_{\boldsymbol{R}\neq\boldsymbol{R}'}\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}\cdots$の部分がややこしいのですが、これは「$\cdots$」の部分を仮に$f(\boldsymbol{R}',\boldsymbol{R})$と置いて考えてみると、

$$
\begin{align*}

\sum_{\boldsymbol{R}\neq\boldsymbol{R}'}\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}f(\boldsymbol{R}',\boldsymbol{R})

&=
\sum_{\boldsymbol{R}}
\sum_{i\in I}
\sum_{\boldsymbol{R}'(\neq\boldsymbol{R})}
\delta_{\boldsymbol{R}',\boldsymbol{R}+\boldsymbol{N}_i}
f(\boldsymbol{R}',\boldsymbol{R})\\

&=
\sum_{\boldsymbol{R}}
\sum_{i\in I}
\sum_{\boldsymbol{R}'(\neq\boldsymbol{R})}
\delta_{\boldsymbol{R}',\boldsymbol{R}+\boldsymbol{N}_i}
f(\boldsymbol{R}',\boldsymbol{R})
\\

&=
\sum_{\boldsymbol{R}}
\sum_{i\in I}

f(\boldsymbol{R}+\boldsymbol{N}_i,\boldsymbol{R})
\\
\end{align*}
$$

と、$\boldsymbol{R}'$を$\boldsymbol{R}+\boldsymbol{N}_i$に変えて、$i\in I$について和を取ればいいことが分かります。

これを踏まえて改めて$\mathcal{H}$を書くと、


$$
\begin{align*}
\mathcal{H}

&=

\sum_{Q_i}\left\{
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    
        \sum_{\boldsymbol{R}}
        
    \varepsilon_{m}^{\rm a}
    \hat{a}_{m,\boldsymbol{R},\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}

    +
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    \sum_{m'\in Q_i}

    \sum_{\boldsymbol{R}\neq\boldsymbol{R}'}
-t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m},\boldsymbol{R})}
\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}
    
    \hat{a}_{m',\boldsymbol{R}',\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}
\right\}
\\

&=

\sum_{Q_i}\left\{
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    
        \sum_{\boldsymbol{R}}
        
    \varepsilon_{m}^{\rm a}
    \hat{a}_{m,\boldsymbol{R},\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}

    +
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    \sum_{m'\in Q_i}

    \sum_{\boldsymbol{R}}
    \sum_{i\in I}
  -t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m},\boldsymbol{R})}
    
    \hat{a}_{m',\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}
\right\}
\\


\end{align*}
$$

ここで例えば、多体シュレーディンガー方程式

$$
\mathcal{H}\Phi = E\Phi
$$

で定まる$N$電子系の固有エネルギー$E$が、大体

$$

\sum_{\boldsymbol{k}\in BZ}\varepsilon_{n_a,\boldsymbol{k}}

\leq

E 

\leq

\sum_{\boldsymbol{k}\in BZ}\varepsilon_{n_b,\boldsymbol{k}}
$$

位の幅であるような状態を考えたいとしましょう。
これは固体内の電子が、特定の$n_a$で指定されるバンドから$n_b$で指定されるバンドの間にとどまっているような状態を想定しています。例えばよくあるのが、価電子が作るバンドだけに電子が詰まっていると考える基底状態を考える場合などです。

これは（多体のハミルトニアンが作用する）スレーター行列式が、原子準位が大体$\varepsilon_{n_a,\boldsymbol{k}}$から$\varepsilon_{n_b,\boldsymbol{k}}$くらいの大きさである場合に対応します。これを適当な幅を$\Delta\varepsilon$として

$$
\varepsilon_{n_a}+\Delta\varepsilon \leq \varepsilon_{m_{q_1}}^{\rm a}, \varepsilon_{m_{q_2}}^{\rm a},
\cdots
,\varepsilon_{m_{q_{|Q|}}}^{\rm a}
\leq
\varepsilon_{n_b} + \Delta\varepsilon
$$

と書いてもそこまで間違っていないと思われます（多分）。

そこで上記を満たすラベルの集合を

$$
Q^{ab}=\{
    m_{q_1^{ab}},m_{q_2^{ab}},\cdots ,m_{q_{|Q^{ab}|}^{ab}}\}

$$

として、第二量子化のハミルトニアンはこの集合に属するラベルのみ和を取る形で

$$
\begin{align*}
\mathcal{H}

&=
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    \sum_{m'\in Q_i}
    
    \sum_{\boldsymbol{R}\boldsymbol{R}'}
    \left\{
    
    \varepsilon_{m}^{\rm a}\delta_{{m}m'}
    \delta_{\boldsymbol{R}'\boldsymbol{R}}

     -t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m},\boldsymbol{R})}
\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}
    \right\}
    \hat{a}_{m',\boldsymbol{R}',\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}\\

&=
\sum_{\gamma=\uparrow,\downarrow}
\left\{
\sum_{m\in Q_i}

\sum_{\boldsymbol{R}}
   
    
    \varepsilon_{m}^{\rm a} 
    \hat{a}_{m,\boldsymbol{R},\gamma}^\dagger
    \hat{a}_{m,\boldsymbol{R},\gamma}

 +
    
    \sum_{m\in Q_i}
    \sum_{m'\in Q_i}
    \sum_{\boldsymbol{R}}
    \sum_{i\in I}
    
    


     -t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m},\boldsymbol{R})}
    \hat{a}_{m',\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}
    \right\}
    \\


\end{align*}
$$

となります。さらに、今このグループの原子準位は大体同じくらいと仮定したので、これを定数$\varepsilon_{ab}^{\rm a}$で置き換えると任意のスレーター行列式に対して

$$
\varepsilon_{ab}^{\rm a}
\sum_{\gamma=\uparrow,\downarrow}
\sum_{m\in Q_i}
\sum_{\boldsymbol{R}}
   
    \hat{a}_{m,\boldsymbol{R},\gamma}^\dagger
    \hat{a}_{m,\boldsymbol{R},\gamma}
    |\cdots|
    =
\varepsilon_{ab}^{\rm a}
|\cdots|
$$

なので第一項を定数とできて、これもエネルギー基準をシフトして消せば、最終的に本稿のゴールであるTight-bindingモデルの第二量子化表示

$$
\begin{align*}
\mathcal{H}

&\simeq
   
\sum_{\gamma=\uparrow,\downarrow}

    
    \sum_{m\in Q_i}
    \sum_{m'\in Q_i}
    \sum_{\boldsymbol{R}}
    \sum_{i\in I}

     -t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m},\boldsymbol{R})}
    \hat{a}_{m',\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}
    \\


\end{align*}
$$

が得られました！

と言ってもまだなんか記号が沢山ありますね。


# 具体例

さて、ようやく第二量子化表示のTight-bindingモデルハミルトニアンが得られました。
ただ、上記の通りまだ何やら記号がたくさんあるので、さらに具体的な例を考えていきます。


## 1つの原子軌道からなるTight-bindingハミルトニアンと対角化

まずは1つだけの原子軌道を考える場合について考えてみます。これは例えば原子準位に縮退のない$s$軌道の関数とその格子ベクトルだけの平行移動$\phi_s(\boldsymbol{r}-\boldsymbol{R})$を考え、まず固有関数がそれらだけで展開できると仮定し、かつ多体の波動関数は$\phi_s(\boldsymbol{r}-\boldsymbol{R})$だけで構成されたスレーター行列式で展開されると仮定することに対応しています。

### 第二量子化表示のハミルトニアン

このように考えた際のハミルトニアンを書いてみると、$s$軌道だけを考えればよいので$\sum_{m\in Q_i}$部分が不要になり

$$
\begin{align*}
\mathcal{H}^s

&=
   
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\boldsymbol{R}}
    \sum_{i\in I}

     -t_{(s,\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({s},\boldsymbol{R})}
    \hat{a}_{s,\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{s,\boldsymbol{R},\gamma}
    \\


\end{align*}
$$

となります。ここで$s$軌道の原子軌道関数をベクトル$\boldsymbol{R}$平行移動した関数$\phi_s(\boldsymbol{r}-\boldsymbol{R})$をスレーター行列式に付け加える（消す）生成（消滅）演算子を、$\hat{a}_{s,\boldsymbol{R},\gamma}^\dagger$（$\hat{a}_{s,\boldsymbol{R},\gamma}$）としました。


この時、大きな特徴として$s$軌道の原子軌道は球対称かつ実関数なので、実関数局所ポテンシャル$V(\boldsymbol{r})$も球対称と仮定すると、飛び移り積分が

$$
\begin{align*}
-t_{(s,\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow (s,\boldsymbol{R})}
&=

\int
 \phi_s(|\boldsymbol{r}-\boldsymbol{N}_i|)
 
   V(|\boldsymbol{r} - \boldsymbol{N}_i|)
 
 \phi_s(|\boldsymbol{r}|)d\boldsymbol{r}

\end{align*}
$$

### 1次元格子、再隣接の飛び移り積分のみ考える

となります。ここで例えば簡単な例として1次元格子を考え、隣接格子ベクトルとして**再隣接格子**を考えてみます。
[第一量子化のTight-bindingモデル（後編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-2)で扱ったように、1次元系を表現する方法は例えば、基本格子ベクトル$\boldsymbol{a}_i$を、

$$
|\boldsymbol{a}_2| = |\boldsymbol{a}_3| = A \gg |\boldsymbol{a}_1|
$$

と置いてしまうことかと思います。これは実際の1次元系について、それ以外の方向との格子間距離が離れている状況を表しています。こうすると$\boldsymbol{a}_2,\boldsymbol{a}_3$方向の再隣接格子の飛び移り積分が、原子軌道関数が中心から離れた場所では十分小さいとしてゼロとできます：

$$
\begin{align*}
-t_{(s,\boldsymbol{R}+\boldsymbol{a}_i) \leftarrow (s,\boldsymbol{R})}
&=

\int
 \phi_s(|\boldsymbol{r}-\boldsymbol{a}_i|)
 
   V(|\boldsymbol{r} - \boldsymbol{a}_i|)
 
 \phi_s(|\boldsymbol{r}|)d\boldsymbol{r}\\

 &\simeq 0\\

 i = 2,3.
\end{align*}
$$

さらにこの時$\boldsymbol{a}_1$方向の再隣接は

$$
\boldsymbol{N}=\pm \boldsymbol{a}_1
$$

ですが、$s$軌道関数が球対称実関数、局所ポテンシャル$V(\boldsymbol{r})$も球対称実関数と置いているので、どちらの方向の飛び移り積分も等しくかつ実数の定数、

$$
\begin{align*}
-t_{(s,\boldsymbol{R}\pm\boldsymbol{a}_1) \leftarrow (s,\boldsymbol{R})}
&=

\int
 \phi_s(|\boldsymbol{r}\pm\boldsymbol{a}_1|)
 
   V(|\boldsymbol{r} \pm \boldsymbol{a}_1|)
 
 \phi_s(|\boldsymbol{r}|)d\boldsymbol{r}\\

 &\equiv-t_s^{\rm NN} = -(t_s^{\rm NN})^*\\

 
\end{align*}
$$

となります。これらを踏まえて第二量子化表示のハミルトニアンは


$$
\begin{align*}
\mathcal{H}^{s,1D}

&=
   
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\boldsymbol{R}}
    \sum_{i\in I}

     -t_{(s,\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({s},\boldsymbol{R})}
    \hat{a}_{s,\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{s,\boldsymbol{R},\gamma}
    \\
&\simeq

\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\boldsymbol{R}}
    \sum_{\boldsymbol{N}_i=\pm \boldsymbol{a}_1}

     -t_{s}^{\rm NN}
    \hat{a}_{s,\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{s,\boldsymbol{R},\gamma}
    \\


\end{align*}
$$


となります。ここで、このハミルトニアンはスレーター行列式の中の1電子状態$\phi_{s,\boldsymbol{R}}$のラベル$\boldsymbol{R}$を、$\boldsymbol{a}_1$方向にしか変えないため、スレーター行列式を1次元上に並んだ原子軌道関数$\phi_s(\boldsymbol{r}), \phi_s(\boldsymbol{r}\pm \boldsymbol{a}_1),\phi_s(\boldsymbol{r}\pm 2\boldsymbol{a}_1),\cdots ,\phi_s(\boldsymbol{r}\pm i\boldsymbol{a}_1)\cdots$のみを考えることにして^[この辺の説明がまだ自分でも混乱している気がするのでいつか修正する]、関数

$$
\phi_s(\boldsymbol{r}\pm i\boldsymbol{a}_1)
$$

（とスピン関数の積）の生成消滅演算子を

$$
    \hat{a}_{s,i,\gamma}^\dagger,\hat{a}_{s,i,\gamma}
$$

と書くことにすると、1次元的なハミルトニアン

$$
\begin{align*}
\mathcal{H}^{s,1D}

&=
   -t_{s}^{\rm NN}
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{i}
    \sum_{\delta=\pm 1}

  
    \hat{a}_{s,i+\delta,\gamma}^\dagger\hat{a}_{s,i,\gamma}
    \\


\end{align*}
$$

が得られます。ここで、$\sum_{i}\sum_{\delta=\pm 1}$ですが、全ての$i$とその隣接項$i\pm 1$についての和を取るので、これを隣接した$i,j$のみ和を取るという意味で、$i<j$として

$$
-t_{s}^{\rm NN}
    \sum_{<i,j>}
    \left(\hat{a}_{s,i,\gamma}^\dagger\hat{a}_{s,j,\gamma}
+
\hat{a}_{s,j,\gamma}^\dagger\hat{a}_{s,i,\gamma}
\right)
$$

と書くこともあります。
余談ですが、この第2項をエルミート共役の意味で"+h.c."と書いていることが多いですが、エルミート共役だと$\hat{a}_{s,j,\gamma}\hat{a}_{s,i,\gamma}^\dagger = -\hat{a}_{s,i,\gamma}^\dagger\hat{a}_{s,j,\gamma}$になってしまうのでつじつまが合わない気がするのですが、正直この辺はあまりよくわかっていません。




### ハミルトニアンの対角化

最後に、このハミルトニアンのシュレーディンガー方程式を解いてみます。つまり固有状態を求めてみます。

それにはハミルトニアンを自明な固有状態を持つ形$\sum_\alpha\varepsilon_\alpha\hat{a}_\alpha^\dagger\hat{a}_\alpha$にすれば良かったのでした。この形を「対角化された」等と呼んだりします。

現状はまだ対角化されていません。


ここで基底変換

$$
\phi_s(\boldsymbol{r}-j\boldsymbol{a}_1)
=
\frac{1}{\sqrt{N}}
\sum_{k_1}
e^{-ik_1a\times j}
\varphi_{s,k_1}(\boldsymbol{r})
$$


を考えると、変換後の$\varphi_{s,k_1}(\boldsymbol{r})$とスピン関数の積の生成消滅演算子を$\hat{b}_{s,k_1,\gamma}^\dagger,\hat{b}_{s,k_1,\gamma}$とすれば、生成消滅演算子の変換は前後でスピン状態は変わらないので

$$
\hat{a}^\dagger_{s,j + \delta,\gamma}
=
\frac{1}{\sqrt{N}}
\sum_{k_1} e^{-i(k_1a\times (j + \delta))}
b_{s,k_1,\gamma}^\dagger,\\

\hat{a}_{s,j,\gamma} = 
\frac{1}{\sqrt{N}}
\sum_{k_1} 
e^{i(k_1a\times j)}
b_{s,k_1,\gamma}
$$

となります。


これをハミルトニアンに代入すると、


$$
\begin{align*}
\mathcal{H}^{s,1D}

&=
   -t_{s}^{\rm NN}
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{j}
    \sum_{\delta=\pm 1}

  \frac{1}{N}\sum_{k_1,k_1'}
  e^{i((k_1-k_1')a\times j)}
    e^{-ik_1'a\delta}
  b_{s,k_1',\gamma}^\dagger
  b_{s,k_1,\gamma}
    \\

&=
   -t_{s}^{\rm NN} \frac{1}{N}
\sum_{\gamma=\uparrow,\downarrow}
    \sum_{\delta=\pm 1}

 \sum_{k_1,k_1'}
  N\delta_{k_1,k_1'}
  e^{-ik_1'a\delta}
  b_{s,k_1',\gamma}^\dagger
  b_{s,k_1,\gamma}
    \\

&=
   -t_{s}^{\rm NN}
\sum_{\gamma=\uparrow,\downarrow}
 \sum_{k_1}
 \sum_{\delta=\pm 1}
  e^{-ik_1a\delta}
  b_{s,k_1',\gamma}^\dagger
  b_{s,k_1,\gamma}
    \\

&=

\sum_{\gamma=\uparrow,\downarrow}
 \sum_{k_1}
 -2t_{s}^{\rm NN}\cos(k_1a)
  b_{s,k_1',\gamma}^\dagger
  b_{s,k_1,\gamma}
    \\
\end{align*}
$$

と、対角化された形になります。これはまさしく、第一量子化で考えた場合と同じく、多電子状態の固有関数はエネルギーバンドが

$$
\varepsilon_s(k) = -2t_{s}^{\rm NN}\cos(ka)
$$

で表される状態（のスピン上向き状態または下向き状態をスレーター行列式に詰めていったもの）となっていることがわかります。



# おわりに

かなり長くなってしまいましたが、以上で第二量子化表示のTight-bindingモデルについての章を終わります。

本章を通して、まずは冒頭に述べた疑問に対して


- サイト$i,j$に電子を「生成する」「消滅させる」って何？そんなことできるの？
  - →**スレーター行列式に付け加える・消す演算子のことでした**
- あとそもそもサイト$i,j$も何？どういう理屈で空間が離散化されてんの？（実はされてません）
- 場の演算子ってやつ？でもそっちは連続変数の位置座標$\boldsymbol{r}$に電子を生成するって書いてあったけど？それを離散化したやつなの？ 
  - →**空間が離散化されているわけではなく、離散的な格子点$\boldsymbol{R}$をラベルとする原子軌道関数$\phi_{\boldsymbol{R}}(\boldsymbol{r})=\phi(\boldsymbol{r}-\boldsymbol{R})$のラベル$\boldsymbol{R}$を指定する添え字なのでした**
- それでその係数の飛び移り積分ってのも何？
  - 隣接格子への飛び移り積分を、隣接格子への相対ベクトルを$\boldsymbol{N}$と置いてのように定義できる積分で、原子軌道関数に固体のハミルトニアンを作用させた場合の展開係数なのでした。

$$
\begin{align*}
-t_{(m,\boldsymbol{R}+\boldsymbol{N}) \leftarrow (n,\boldsymbol{R})}
&\equiv

\int
 \phi_m^*(\boldsymbol{r}-\boldsymbol{N})
 
   V(\boldsymbol{r} - \boldsymbol{N})
 
 \phi_n(\boldsymbol{r})d\boldsymbol{r}\\

 &=

 \int
 \phi_m^*(\boldsymbol{r}-(\boldsymbol{R}+\boldsymbol{N}))
 
   V(\boldsymbol{r} - (\boldsymbol{R}+ \boldsymbol{N}))
 
 \phi_n(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}

\end{align*}
$$
- ていうか全部何？
  - →本章で学んだことを以下にまとめておきます

## 本章のまとめ

また本章の概要をまとめておくと、まず第二量子化表示の概要として、

- (1) 多体電子系のハミルトニアンを設定して1体演算子の総和で表し、
- (2) スレーター行列式を展開するための何かしらの完全正規直交関数系を決めて
- (3) 完全正規直交関数系でハミルトニアンを挟んで積分したハミルトニアン行列を求める

ことで、多体のハミルトニアンの任意のスレーター行列式に対する作用として、

$$
\mathcal{H}\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|=\sum_i \hat{H}_i
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
    
=    
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

と生成消滅演算子を用いて表現でき、これを第二量子化表示というのでした。特にTight-bindingモデルにおいては、上記3ステップは以下のようになります。

### (1) 多体電子系のハミルトニアンを設定する：
:::message
孤立原子のポテンシャルが、固体内の格子点に周期的に並んでいると考え、かつ電子間のCoulomb相互作用を無視して

$$
\mathcal{H}

\simeq

\sum_i
  \left(
\frac{-\hbar^2}{2m}\nabla_i{}^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}_i-\boldsymbol{R}'')
\right) 

\equiv
\sum_i\hat{H}_i^{\rm c}
$$

と一体近似をしたハミルトニアンを考える。

:::
### (2) スレーター行列式を展開するための何かしらの完全正規直交関数系を決める：

:::message
固体のハミルトニアンを$\hat{H}^{\rm c}=\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})$と1体近似で近似した際、周期的に並んでいるとしたポテンシャル$V(\boldsymbol{r})$が「孤立して」ある場合のハミルトニアン（孤立原子のハミルトニアン）$\hat{H}^{\rm a}=-(\hbar^2/2m)\nabla^2 + V(\boldsymbol{r})$の固有関数：

$$
\hat{H}^{\rm a} \phi_m(\boldsymbol{r}) = \varepsilon_m^{\rm a}\phi_m(\boldsymbol{r}). 
$$

およびその格子ベクトル分平行移動した関数系

$$
\phi_m(\boldsymbol{r}-\boldsymbol{R})
$$

なお、LCAO近似・Tight-binding近似より、この原子軌道関数のうちエネルギー的に近い原子準位を持つ原子軌道のラベルの集合

$$
Q = \{m_{q_1},m_{q_2},\cdots m_{q_{|Q|}}\}
$$

を考え、$Q$に属する原子軌道の原子準位を$\varepsilon_{m_{q_1}}^{\rm a}, \varepsilon_{m_{q_2}}^{\rm a},\cdots \varepsilon_{m_{q_{|Q|}}}^{\rm a}$とし、これら原子準位と近い固有値をもつ固体のハミルトニアンの固有状態$\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r})$：$\hat{H}^{\rm c}\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r}) = \varepsilon_{n_{p_i},\boldsymbol{k}}\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r})$、$\varepsilon_{n_q,\boldsymbol{k}} \simeq \varepsilon_{{q_1}}^{\rm a}, \cdots$とする。

すると$|Q|$個の原子軌道関数から作られたBloch和の重ね合わせから、$|Q|$個の固体の固有関数を以下のように**近似的に**展開できる：

$$
\varphi_{n_{p_i},\boldsymbol{k}}(\boldsymbol{r}) \simeq \sum_{m \in Q}b_m^{n_{p_i}}\sum_{\boldsymbol{R}}e^{i\boldsymbol{k}\cdot\boldsymbol{R}} \phi_m(\boldsymbol{r}-\boldsymbol{R}),\\
i = 1,2,\cdots |Q|.
$$

また、この展開によって得られる固有関数のラベルの集合を

$$
P = \{n_{p_1}, n_{p_2}, \cdots \} 
$$

と置き、$|P|=|Q|$を満たす。


:::
  
### (3) 完全正規直交関数系でハミルトニアンを挟んで積分したハミルトニアン行列を求める

:::message
上記の完全正規直交関数系と、固有関数の近似的な展開の関係を用いて、

$$
\begin{align*}
&\braket{m'\boldsymbol{R}'|\hat{H}^{\rm c}|{m_{q_j}}\boldsymbol{R}}\\
&=\int \phi_{m'}^*(\boldsymbol{r}-\boldsymbol{R}')\hat{H}\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
\\
&=
\int\phi^*_{m'}(\boldsymbol{r} - \boldsymbol{R}')\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}-\boldsymbol{R}'')
\right)\phi_{m_{q_j}}(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
\\

&\simeq

\left\{
\left(
\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}

+\Delta\varepsilon_{m'{m_{q_j}}}
\right)\delta_{\boldsymbol{R}'\boldsymbol{R}}


-t_{(m',\boldsymbol{R}') \leftarrow ({m_{q_j}},\boldsymbol{R})}(1-\delta_{\boldsymbol{R}'\boldsymbol{R}})
\right\}\delta_{m'\in Q}'\\

&\simeq
\left\{
\left(
\varepsilon_{m_{q_j}}^{\rm a}\delta_{{m_{q_j}}m'}

+\Delta\varepsilon_{m'{m_{q_j}}}
\right)\delta_{\boldsymbol{R}'\boldsymbol{R}}


-t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m_{q_j}},\boldsymbol{R})}
\sum_{i\in I}\delta_{\boldsymbol{R}'-\boldsymbol{R},\boldsymbol{N}_i}
\right\}\delta_{m'\in Q}'
\end{align*}
$$

となる。
ここで、飛び移り積分を

$$
-t_{(m,\boldsymbol{R}') \leftarrow (n,\boldsymbol{R})}
\equiv

\int
 \phi_m^*(\boldsymbol{r}-\boldsymbol{R}')
 
   V(\boldsymbol{r} - \boldsymbol{R}')
 
 \phi_n(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}

$$

定義した。最後の近似は隣接格子のみの飛び移り積分を考える近似をすると、隣接格子への飛び移り積分を、隣接格子への相対ベクトルを$\boldsymbol{N}$と置いて

$$
\begin{align*}
-t_{(m,\boldsymbol{R}+\boldsymbol{N}) \leftarrow (n,\boldsymbol{R})}
&\equiv

\int
 \phi_m^*(\boldsymbol{r}-\boldsymbol{N})
 
   V(\boldsymbol{r} - \boldsymbol{N})
 
 \phi_n(\boldsymbol{r})d\boldsymbol{r}\\

 &=

 \int
 \phi_m^*(\boldsymbol{r}-(\boldsymbol{R}+\boldsymbol{N}))
 
   V(\boldsymbol{r} - (\boldsymbol{R}+ \boldsymbol{N}))
 
 \phi_n(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}

\end{align*}
$$

と書け、どちらも前後の格子点の座標の差のみに依存する。
:::

### Tight-bindingモデルの第二量子化表示

:::message

完全正規直交関数系に対応して、関数$\phi_{m,\boldsymbol{R}}(\boldsymbol{r})=\phi_{m}(\boldsymbol{r}-\boldsymbol{R})$とスピン関数$\gamma(\sigma)=\alpha(\sigma), \beta(\sigma)$の積

$$
\phi_{m,\boldsymbol{R}}(\boldsymbol{r})\gamma(\sigma)
\equiv
\phi_{m,\boldsymbol{R},\gamma}(\boldsymbol{r})
$$

をスレーター行列式に「付け加える」（既にあれば$0$を返す）生成演算子を$\hat{a}_{m,\boldsymbol{R},\gamma}^\dagger$：

$$
\hat{a}_{m,\boldsymbol{R},\gamma}^\dagger|\cdots|=|\phi_{m\boldsymbol{R},\gamma}\cdots|
$$

スレーター行列式の先頭から$\phi_{m,\boldsymbol{R},\gamma}(\boldsymbol{r})$を消す（無ければ$0$を返す、先頭ではない場所にあれば先頭まで移動してから消す）消滅演算子を$\hat{a}_{m,\boldsymbol{R},\gamma}$：

$$

\hat{a}_{m,\boldsymbol{R},\gamma}|\phi_{m,\boldsymbol{R},\gamma}\cdots|=|\cdots|

$$

をを定義する。


結晶場積分をラベル$m,m'$によらない定数と近似し、基準となるエネルギーを結晶場積分の値の分だけずらして式から消去して、かつ飛び移り積分を隣接格子だけ取ることにすると、Tight-bindingモデルの第二量子化表示は、あるラベル$m$が属する集合を$Q(\ni m)$と書くことにして、さらに全ての$m$の集合を$Q_1, Q_2\cdots Q_i \cdots$と分割することにすると、（こんな書き方が正しいかどうかはわからないのですが）


$$
\begin{align*}
\mathcal{H}



&=

\sum_{Q_i}\left\{
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    
        \sum_{\boldsymbol{R}}
        
    \varepsilon_{m}^{\rm a}
    \hat{a}_{m,\boldsymbol{R},\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}

    +
    \sum_{\gamma=\uparrow,\downarrow}
    \sum_{m\in Q_i}
    \sum_{m'\in Q_i}

    \sum_{\boldsymbol{R}}
    \sum_{i\in I}
  -t_{(m',\boldsymbol{R}+\boldsymbol{N}_i) \leftarrow ({m},\boldsymbol{R})}
    
    \hat{a}_{m',\boldsymbol{R}+\boldsymbol{N}_i,\gamma}^\dagger\hat{a}_{m,\boldsymbol{R},\gamma}
\right\}
\\


\end{align*}
$$

と、**「エネルギーが近いグループ」ごとに分かれた形**になる。
:::


というわけでようやくゴールに到達です。まだ若干ごまかしている部分や、もう少し追記したい部分もあったり、2体演算子についても考えていきたいところですが、そのあたりはぼちぼちまたやることにして、ひとまず一人で学部卒業式でも敢行したいと思います。

