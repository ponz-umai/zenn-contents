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
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }\left( -t_{ij}a^\dagger_{i\sigma }a_{j\sigma }   + h.c.\right)
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

- 多体電子系のハミルトニアンを設定する
- スレーター行列式を展開するための何かしらの完全正規直交関数系を決める

ことだとわかります。この二つが決まれば、後はこれまでに求めた手続きに従うことで第二量子化表示ができそうです。なお、二つ目の「完全正規直交関数系」を「基底」、「基底関数」などと呼び、特定の関数を選ぶことを「基底を決める」などと呼びます。この辺はベクトルっぽい用語が示すように「関数のベクトル表現」とかいう話と密接に関わる話なのですが、本稿では~~ギリ触れなくても説明できそうな気がしたので~~触れていません。が、「基底」という用語は字数も少なく済むので適宜使用します。^[皆さんご存じの内容とは思いますが、「関数のベクトル表現」「Hilbert空間」「Fock空間」などのキーワードで調べたり、適当な量子力学の教科書を読んだりするとちゃんと勉強できると思います。]

早速一つ目から始めていきましょう。

## 多体電子系のハミルトニアン：固体中の多体電子系のハミルトニアンを設定する


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

## スレーター行列式を展開する完全正規直交関数系を決める

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

つまり、「完全正規直交関数系」第一候補は平面波$\{\frac{1}{\sqrt{V}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}\}$です。

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

によって書くことができます。このように固有関数は、組み合わせる波数を代表する波数$\boldsymbol{k}$（この波数を「Bloch波数」と呼びます）と、（無限個の）逆格子ベクトルだけ異なる平面波の組み合わせから得られる（無限個の）固有値に対応したラベル$n$によって指定され、Bloch関数などと呼ばれます。Bloch関数に対応する固有値を、ラベル$n,\boldsymbol{k}$を用いて$\varepsilon_{n,\boldsymbol{k}}$と書き、固有値、固有関数（Bloch関数）は$\hat{H}^{\rm c}\varphi_{n,\boldsymbol{k}} = \varepsilon_{n,\boldsymbol{k}}\varphi_{n,\boldsymbol{k}}$を満たします。固有値のラベル$n$は、固有値（固有エネルギー）が小さい順に番号を振ることにしておきます。

#### 完全正規直交基底その2：Bloch関数

さて、前置きが長くなりましたが、上記のようにして構成されるBloch関数は正規直交性$\int\varphi^*_{n',\boldsymbol{k}'}(\boldsymbol{r})\varphi_{n,\boldsymbol{k}}(\boldsymbol{r})d\boldsymbol{r} = \delta_{n',n}\delta_{\boldsymbol{k}',\boldsymbol{k}}$を満たしまた、完全正規直交性を満たす平面波のセット$\{\frac{1}{\sqrt{V}}e^{i\boldsymbol{q}\cdot\boldsymbol{r}}\}$の線形結合、いわば組み換えなわけなので任意の関数を表すことができる完全性も満たします。

:::details 一応証明しておきます
少し抽象的になってしまいますが、平面波のセット$\{e^{i\boldsymbol{q}\cdot\boldsymbol{r}}\}$とBloch関数のセット$\{\varphi_{n,\boldsymbol{k}}\}$は、
Fourier係数が満たす方程式

$$
 (\frac{\hbar^2}{2m}q^2 - \varepsilon)c_{\boldsymbol{q}} + 
    \sum_{\boldsymbol{K}}U_{\boldsymbol{K}}
c_{\boldsymbol{q}-\boldsymbol{K}} = 0
$$

によって結び付けられていますが、

$$
(H^{\boldsymbol{q}})_{ij} \equiv \hbar^2(\boldsymbol{q}-\boldsymbol{K}_i)^2/2m\delta_{ij} + U_{\boldsymbol{K}_i+\boldsymbol{K}_j}
$$

とベクトル

$$
\boldsymbol{c}^{\boldsymbol{q}}=(\cdots , c_{\boldsymbol{q} - \boldsymbol{K}_{i+1}}, c_{\boldsymbol{q} - \boldsymbol{K}_i},\cdots c_{\boldsymbol{q}},\cdots, c_{\boldsymbol{q} + \boldsymbol{K}_i}, c_{\boldsymbol{q} + \boldsymbol{K}_{i+1}},\cdots)^T
$$

を定義した行列表示//こんなんせんでも、c = <e|phi>として変換行列を定義して逆行列を定義すれば良い

$$

$$

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

と書くこともできます。（この式を利用する方が多いかも）

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


先ほど「Wannier関数は実空間上の位置座標$\boldsymbol{r}$の関数であり、かつ実空間上の格子点の座標でラベルされるややこしい見た目をしています。」と書きましたが、ラベル$\boldsymbol{R}$は実空間上の関数の「中心」または「平行移動する分」を意味しているのですね。

そして、このようにして定義されたWannier関数

$$
\begin{align*}
w_{n,\boldsymbol{R}}(\boldsymbol{r})&=

\frac{1}{v_{BZ}}\int_{BZ} \varphi_n(\boldsymbol{r},\boldsymbol{k})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}d\boldsymbol{k}\\
&=
\frac{1}{N}\sum_{\boldsymbol{k}\in BZ} \varphi_{n\boldsymbol{k}}(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
\end{align*}
$$

もまた、完全正規直交関数系であるBloch関数$\{\varphi_{n\boldsymbol{k}}\}$から組み替えられたものであり、完全正規直交性を満たす基底になり得ます。

とはいえそもそものBloch関数の具体的な関数系が指定されなければ、Wannier関数やWannier関数によるハミルトニアンの積分も計算できません。次へ進みましょう。


#### （近似的な）完全正規直交関数：原子軌道関数

さて、前置きが長くなりましたが、ようやく具体的な関数系を導入します。

Wannier関数の「局在性」に注目し、原点に局在するWannier関数を、同じく原点を中心とする完全正規直交関数系である「原子軌道関数」を用いて展開しよう、と考えるわけです。ここで原子軌道関数を以下のように、孤立原子ポテンシャル1つだけがある系の固有関数として定義します。（孤立原子ポテンシャルは固体中のポテンシャルを$\hat{H}^{\rm c}=\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r} - \boldsymbol{R})$とを1体近似で近似した際、周期的に並んでいるとしたポテンシャル$V(\boldsymbol{r})$です。）



すなわち、孤立原子のハミルトニアンを$\hat{H}^{\rm a}=-(\hbar^2/2m)\nabla^2 + V(\boldsymbol{r})$として、

$$
\hat{H}^{\rm a} \phi_m(\boldsymbol{r}) = \varepsilon_m^{\rm a}\phi_m(\boldsymbol{r}). 
$$

ここで固有関数$\phi_m(\boldsymbol{r})$を「原子軌道」や「原子軌道関数」、固有値$\varepsilon_m^{\rm a}$を「原子準位」等と呼ぶことにします。

このように定義した原子軌道関数を用いて、先ほどBloch関数を展開した「格子点$\boldsymbol{R}$に中心を持つ」Wannier関数$w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_{n}(\boldsymbol{r}-\boldsymbol{R})$を、「格子点$\boldsymbol{R}$に中心を持つ」孤立原子の固有関数$\phi_m(\boldsymbol{r}-\boldsymbol{R})$と展開係数$b_m$を用いて、

$$
w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_n(\boldsymbol{r}-\boldsymbol{R}) = \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R}) 
$$

と展開します。ここで、Wannier関数の並進性を考えると任意の$\boldsymbol{R}$ラベルの$w_{n,\boldsymbol{R}}$は、別のラベル$\boldsymbol{R}'$を平行移動しただけの関数ですので、展開係数$b_m^n$は$\boldsymbol{R}$によらずただ原子軌道関数が$\boldsymbol{R}$だけ平行移動されるだけとなります。

また規格化条件$\sum_m|b_m^n|^2 = 1$を満たすものとします。

また、逆変換を

$$
w_{n,\boldsymbol{R}}(\boldsymbol{r}) = w_n(\boldsymbol{r}-\boldsymbol{R}) = \sum_mb_m^n\phi_m(\boldsymbol{r}-\boldsymbol{R}) 
$$

このような変換を考えると、完全正規直交関数系Wannier関数


その代わり、原子軌道間の**重なり積分**を以下のように定義し、



//この辺もうちょいうまい説明が欲しい











これは、前章の第二量子化の考え方に沿って言い直すと、各格子点を中心とする原子軌道関数$\{\phi_m(\boldsymbol{r} - \boldsymbol{R})\}$または$\boldsymbol{R}$もラベルとして書いて$\{\phi_{m,\boldsymbol{R}}(\boldsymbol{r})\}$を軌道関数部分の完全正規直交基底として、この関数系とスピン関数$\gamma=\alpha,\beta$どちらかの積の完全正規直交系のスピン軌道関数を$\{\phi_{m,\boldsymbol{R}, \sigma}(\boldsymbol{r})\}$を並べて構成したスレーター行列式$|\phi_{n,\boldsymbol{R}_i,\sigma_p}, \phi_{l,\boldsymbol{R}_j,\sigma_q},\cdots\phi_{o,\boldsymbol{R}_k,\sigma_r}|$で多電子状態の波動関数を

$$
\Phi(\tau_1, \tau_2, \cdots, \tau_N) = \sum C_{}\left|
    \phi_{n,\boldsymbol{R}_i,\sigma_p}, \phi_{l,\boldsymbol{R}_j,\sigma_q},\cdots
    \phi_{o,\boldsymbol{R}_k,\sigma_r}
    \right|
$$

と展開すると考えることに対応する。

この時前章のレシピを利用して、特に1体のハミルトニアン$\frac{-\hbar^2}{2m}\nabla_i{}^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}_i - \boldsymbol{R})$にスピンを含む項がないため、


$$
\begin{align*}
\mathcal{H} &= \sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\\

&=
\sum_{n,\boldsymbol{R}_i;l,\boldsymbol{R}_j \gamma=\alpha,\beta}
\int \phi_{n,\boldsymbol{R}_i}^*(\boldsymbol{r})\hat{H}\phi_{l,\boldsymbol{R}_j}(\boldsymbol{r})d\boldsymbol{r}
\hat{a}^\dagger_{n,\boldsymbol{R}_i\gamma} \hat{a}_{l,\boldsymbol{R}_j\gamma}\\

&\equiv
\sum_{n,\boldsymbol{R}_i;l,\boldsymbol{R}_j \gamma=\alpha,\beta}
\braket{n,\boldsymbol{R}_i|\hat{H}|l,\boldsymbol{R}_j}
\hat{a}^\dagger_{n,\boldsymbol{R}_i\gamma} \hat{a}_{l,\boldsymbol{R}_j\gamma}\\

\end{align*}
$$

と書ける。そこで残すはハミルトニアン行列要素$\braket{n,\boldsymbol{R}_i|\hat{H}|l,\boldsymbol{R}_j} = \int \phi_{n,\boldsymbol{R}_i}^*(\boldsymbol{r})\hat{H}\phi_{l,\boldsymbol{R}_j}(\boldsymbol{r})d\boldsymbol{r}$の計算となりました。

この積分はどこかで見覚えがあるものです、どこかというと[飛び移り積分（Hopping Integral）の物理的意味・Wannier関数の従う方程式](https://zenn.dev/ponzumai/articles/tight-binding-model-hopping-int)の章で、Tight-binding近似の下では、原子軌道関数に固体のハミルトニアンを作用させた結果は

$$
\begin{align*}
\hat{H}^{\rm c} \phi_m(\boldsymbol{r})
&=
\left(
\frac{-\hbar^2}{2m}\nabla^2 + \sum_{\boldsymbol{R}''}V(\boldsymbol{r}-\boldsymbol{R}'')
\right) \phi_m(\boldsymbol{r}) \\


&\simeq

(\varepsilon_m^{\rm a} + \Delta\varepsilon_{mm})\phi_m(\boldsymbol{r})

+
\sum_{m'} \Delta\varepsilon_{m'm}\phi_{m'}(\boldsymbol{r})\\

&\>\>\>\>+\sum_{\boldsymbol{R}'}\sum_{m'}(-t_{\boldsymbol{R}'}^{m,m'})^*\phi_{m'}(\boldsymbol{r} - \boldsymbol{R}')

\end{align*}
$$

となり、これに格子点$\boldsymbol{R}$


～～～～

というわけでこれまでの問題設定の下多体のハミルトニアンを第二量子化表示すると

# 具体例

## 1つの原子軌道からなるTight-bindingハミルトニアンと対角化

## 縮退のある原子軌道からなるTight-bindingハミルトニアン

