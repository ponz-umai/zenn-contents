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
- それでその係数の飛び移り積分ってのも何？
- あとどういう理屈で空間が離散化されてんの？（実はされてません）
- 場の演算子ってやつ？でもそっちは連続変数の位置座標$\boldsymbol{r}$に電子を生成するって書いてあったけど？それを離散化したやつなの？（違います）

等と浮かんでは消える疑問から目を背けながら、目の前に積みあがる研究タスクを消化するために「反交換関係」とやらを駆使しながら過ごす日々をお過ごしの方が、きっとこのような怪しげな個人ブログに行きついてきているものと思います。（そうですよね？それとも↑のような疑問に頭を悩ませていたのは私だけ？）

そんな日々に終止符を打つため、本章ではこれまでまとめてきた内容をもとに、上記のような第二量子化表示のTight-bindingモデルについて解説していきます。

# 前章の振り返り

初めに今回使う道具として、[前章](https://zenn.dev/ponzumai/articles/tight-binding-model-2nd-q)でまとめた第二量子化表示の基本的な考え方をこちらにまず再掲します。

まず、多体電子系のシュレーディンガー方程式を

$$
\mathcal{H}\Phi(\tau_1,\tau_2,\cdots,\tau_N)
=
E
\Phi(\tau_1,\tau_2,\cdots,\tau_N)
$$

と書いた上で、任意の多体の波動関数$\Phi(\tau_1,\tau_2,\cdots,\tau_N)$をスレーター行列式の線形結合として

$$
\Phi(\tau_1,\tau_2,\cdots,\tau_N)
=
\sum_{\lambda,\mu,\cdots,\xi}C_{\lambda,\mu,\cdots\xi}
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

と展開すると考えると、多体のハミルトニアン$\mathcal{H}$が1電子のハミルトニアン$\hat{H}$の和で

$$
\mathcal{H} = \sum_i\hat{H}_i
$$

と書かれている場合ハミルトニアンの作用は

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

となります。

ここで、以下のように生成演算子$\hat{a}_\delta^\dagger$と消滅演算子$\hat{a}_\kappa$を定義することで：

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

として変換できます。

また、特に1体のハミルトニアンがスピンに関する項を含まない場合、第二量子化表示のハミルトニアンは


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

# 問題設定

さて、早速Tight-bindingモデルの問題設定から始めます。これも基本的にこれまでの章、特に第一量子化のTight-bindingモデル[（前編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-1)[（後編）](https://zenn.dev/ponzumai/articles/tight-binding-model-1st-q-2)で扱ってきた内容の復習です。

まずは上記の章で述べたように、大きな目標としては固体、すなわち周期的な格子ポテンシャルが存在するような空間の中に、多数の電子がいるような物理系における電子の振る舞い、すなわち波動関数を求めたいわけです。

どんなややこしい状況であれ電子の波動関数を求める基本方程式は、考えている系に存在するポテンシャルをもとに構成された多体のハミルトニアン$\mathcal{H}$と、多体の波動関数$\Phi(\tau_1, \tau_2,\cdots\tau_N)$によるシュレーディンガー方程式

$$
\mathcal{H}\Phi(\tau_1, \tau_2,\cdots\tau_N) = E\Phi(\tau_1, \tau_2,\cdots\tau_N)
$$

です。しかしながら、ありのまますべて正確に書こうとすると、固体内の相互作用はこんな感じで


![](/images/tb/many-ele-atm.png)

多数の原子核の相互作用を受けながら、多数の電子がそれぞれに相互作用をしあっている状態です。わけわからんですね。ハミルトニアンもそれに対応して、

（原子核（や分子）からの相互作用ポテンシャル）＋（原子核同士の相互作用ポテンシャル）＋（電子間の相互作用ポテンシャル）＋（その他スピンに関する相互作用やらなんやかんや）

のように到底解けない形になってしまうわけです。

## ポテンシャルの近似

そこで、まず第一の近似として、

- 原子核（格子点）は静止しているものと近似する（Born–Oppenheimer近似）
- 内殻電子と価電子間の相互作用は全て価電子が感じる1体ポテンシャルとして近似し、結晶の周期と同じ周期で局所的な孤立原子のポテンシャルとして扱う
- 価電子間のCoulomb相互作用はひとまず無視する

のようにして、**お互いに相互作用しない価電子が、周期的に並んだ（静止した）局所ポテンシャル（孤立原子ポテンシャル）を受けながら運動する**モデルを考えることにするわけです。これを1体近似などと呼びます。

これにより、孤立原子ポテンシャルを$V(\boldsymbol{r})$、結晶の周期を表す基本格子ベクトルを$\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3$と書き、各格子点の位置ベクトルを$\boldsymbol{R} = n_1\boldsymbol{a}_1 + n_2\boldsymbol{a}_2 + n_3\boldsymbol{a}_3, n_i = 0,\pm 1, \pm 2\cdots$と書いて、多体のハミルトニアンを

$$
\mathcal{H} \simeq \sum_i\left\{

\frac{-\hbar^2}{2m}\nabla_i{}^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}_i - \boldsymbol{R})

    \right\}
$$

と近似します。これでハミルトニアンを前章で整理した1体演算子の総和の形に書くことができました。


## Blochの定理・周期ポテンシャル下の固有関数とWannier関数

1体演算子の総和の時はスレーター行列式の中の固有関数が1体ハミルトニアンの固有関数になる？みたいなことを書く？

～～～～～

ここで、[周期ポテンシャル内の電子状態の章](https://zenn.dev/ponzumai/articles/tight-binding-model-electrons-in-solids)上記のような周期ポテンシャルの下で、電子が取る固有関数はBloch関数として表され、さらに[Wannier関数の章](https://zenn.dev/ponzumai/articles/tight-binding-model-wannier-func)で見たように、Bloch関数は局在したWannier関数$w_{\boldsymbol{R}}(\boldsymbol{r}) = w(\boldsymbol{r} - \boldsymbol{R})$の線形結合

$$
\varphi_{\boldsymbol{k}}(\boldsymbol{r}) =\sum_{\boldsymbol{R}}w_{\boldsymbol{R}}(\boldsymbol{r})e^{i\boldsymbol{k}\cdot\boldsymbol{R}},\\
w_{\boldsymbol{R}}(\boldsymbol{r})
=
\frac{1}{N}\sum_{\boldsymbol{k}} \varphi_{\boldsymbol{k}}(\boldsymbol{r})e^{-i\boldsymbol{k}\cdot\boldsymbol{R}}
$$


により展開できました。

またWannier関数は
- 並進性：$w_{\boldsymbol{R}}(\boldsymbol{r}) = w(\boldsymbol{r}-\boldsymbol{R})$
- 正規直交性：$\int_V w^*(\boldsymbol{r}-\boldsymbol{R}')w(\boldsymbol{r} - \boldsymbol{R})dr =\delta_{\boldsymbol{R},\boldsymbol{R}'}$
- 局在性：ラベル＝格子点の座標$\boldsymbol{R}$を中心として局在した関数

を満たすのでした。なお、ここではバンド指標$n$は指定せずに書いています。


## LCAO近似

ここでさらに、孤立原子ポテンシャルが周期的に並んだ系を考えていることから、**電子の波動関数もこの孤立原子の固有関数「原子軌道関数」の重ね合わせで書けるだろう**と考える**LCAO(Linear Combination of Atomic Orbitals)近似**を導入します。

ここで、孤立原子の固有関数$\phi_m(\boldsymbol{r})$を以下のように定義します：

孤立原子のハミルトニアンを$\hat{H}^{\rm a}=-(\hbar^2/2m)\nabla^2 + V(\boldsymbol{r})$として、

$$
\hat{H}^{\rm a} \phi_m(\boldsymbol{r}) = \varepsilon_m^{\rm a}\phi_m(\boldsymbol{r}). 
$$

ここで固有関数$\phi_m(\boldsymbol{r})$を「原子軌道」や「原子軌道関数」、固有値$\varepsilon_m^{\rm a}$を「原子準位」等と呼ぶことにします。

このように定義した原子軌道関数を用いて、先ほどBloch関数を展開したWannier関数を、孤立原子（又は分子等）の波動関数$\phi_m(\boldsymbol{r})$と展開係数$b_m$を用いて、

$$
w_{\boldsymbol{R}}(\boldsymbol{r}) = w(\boldsymbol{r}-\boldsymbol{R}) = \sum_mb_m\phi_m(\boldsymbol{r}-\boldsymbol{R}) 
$$

と展開します。$b_m$は規格化条件$\sum_m|b_m|^2 = 1$を満たすものとします。

ここであるバンドに対応する（Bloch関数を展開した）Wannier関数を、少数の原子軌道で展開できる、と仮定（近似）するのがLCAO近似ですが、ここではまだその近似は取り入れずにおきます。

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

