---
title: "第二量子化表示/生成・消滅演算子の方法"
emoji: "🚥"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# はじめに

前章までで（ポエムを挟みつつ）ようやく第一量子化でのTight-bindingモデルの解説が終わりました。

最終ゴールであるTight-bindingモデルの第二量子化表示へ向け、本章では多体電子系の波動関数、そしてハミルトニアンを扱うための手法である、第二量子化表示について解説していきます。

第二量子化表示は、（皆さんご存じの通り）

$$
\mathcal{H}


=
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa + \cdots
$$

のように、多体電子系の波動関数、そしてハミルトニアンを「生成・消滅演算子」と呼ばれる演算子を用いて書き表す方法です。

第二量子化表示を見て、「電子をそんな勝手に作ったり消したりしていいの？」とか「この係数って何？」とか思い悩んだり、「場の演算子」とかいうワードにビビったりなどして、何かすごい難しいものだという印象を受ける第二量子化ですが、本章では厳密性なことは諦めつつ、何をどう考えればそういう風に書き直せるのか、可能な限り丁寧に書いていこうと思います。

また、最近の本やWEBにある解説では場の演算子から出発するものが多いですが、個人的に好みである（おそらく古めの本に多い）導入方法、すなわち、「第一量子化の」波動関数やハミルトニアン、シュレーディンガー方程式を第二量子化表示へ書き換えていく導入方法で、第二量子化を解説していきます。（場の演算子については触れません）

ただ、例によって変なことを書いているかもしれませんので、ちゃんと勉強したい方は世にあるまじめな解説^[例えばWEBで読める記事だと[田崎先生](https://www.gakushuin.ac.jp/~881791/pdf/2ndQNoteJ.pdf)や[加藤先生](https://kato.issp.u-tokyo.ac.jp/kato/lecture_ISSP/second_quantization_note_v5.pdf)や[山﨑先生](https://home.hiroshima-u.ac.jp/kyam/pages/results/monograph/Ref29_Fock.pdf)、[永井先生による動画での解説](https://www.youtube.com/watch?v=Gl6yVPU12Wc)なんかもあります^[物理界隈だと、「年上であっても研究者のことは「さん」付けで呼ぶのが分かっている奴であり、「先生」と呼ぶのはダサい」みたいな空気があるが、とはいえ上記先生方と面識があるわけでもないのでここは「先生」としておく]。教科書だと、私が参考にしたのが小出先生の[『量子力学(II)』](https://www.amazon.co.jp/dp/4785321431/)（最近新装版が発刊されたみたいですね）です。またWEBの記事でも紹介した加藤先生の[『一歩進んだ理解を目指す　物性物理学講義』](https://www.amazon.co.jp/dp/4781915337)では、場の演算子から出発する考え方や、第二量子化の簡単な応用も丁寧に解説されていました。]を是非参照してください。


# 多体電子系における1体演算子の第二量子化表示

早速始めていきましょう。第二量子化について勉強し始めると、やれ場の演算子だの、生成・消滅演算子の定義だの細かい話が出てきてそこで一旦教科書・ブラウザを閉じることが多いかと思います。ということを避けるために、ここでは個人的に一番面白いというか、ああ、そういうことねと思ったところから話を始めようと思います。それがスレーター行列式への演算子の作用についてです。

## イントロダクション：第二量子化の雰囲気

これまでに扱ってきたように、多体のハミルトニアン$\mathcal{H}$が1電子のハミルトニアン$\hat{H}$の和で

$$
\mathcal{H} = \sum_i\hat{H}_i
$$

と書かれているとしましょう。ここで$\hat{H}_i$は、例えば固体中の電子のハミルトニアンのように、$i$番目の電子の座標を$\boldsymbol{r}_i$、これへの微分演算子を$\nabla_i$として、$\hat{H}_i = -(\hbar^2/2m)\nabla_i^2 + \sum_{\boldsymbol{R}}V(\boldsymbol{r}_i -\boldsymbol{R})$や、あるいは本稿では扱っていませんが$i$番目のスピンに対する作用（ポテンシャル）を含んでいても良く、とにかく1電子だけにかかわる演算子であればなんでも良いです。

この演算子を$N$電子の波動関数に対応するスレーター行列式$\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|$に作用させると、


$$
\begin{align*}
&\mathcal{H}\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|\\



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
    \right|

\end{align*}
$$

と、スレーター行列式の各1電子の波動関数のどれか一つに対して、1体ハミルトニアンが作用したスレーター行列式の総和となります。（証明は後でやります）


ここで軌道部分の波動関数$\varphi_a, \varphi_b\cdots$は正規直交性を満たすものとし、また完全系を成すと仮定します。（スピン関数もスピン空間内で完全系を成すので、その積のスピン軌道関数も完全系をなすと仮定していることになります）

すると1電子波動関数$\varphi_\lambda, \varphi_\mu,\cdots\varphi_\xi$は完全系であると仮定したので、それぞれに1体のハミルトニアンが作用した後の関数を$\hat{H}\varphi_\kappa(\tau) = \sum_{\delta}c_{\delta}^{\kappa}\varphi_{\delta}(\tau)$と展開でき、また正規直交性より、展開係数$c_{\delta}^{\kappa} = \int \varphi_{\delta}^*(\tau)\hat{H}\varphi_{\kappa}(\tau)d\tau$となるので、$\int \varphi_{\delta}^*(\tau)\hat{H}\varphi_{\kappa}(\tau)d\tau\equiv H_{\delta\kappa}$と書くと、（$\tau$は[電子のスピンを考慮した多電子系の波動関数](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)で書いたように、座標の変数$\boldsymbol{r}$とスピン変数$\sigma$を合わせた変数です。細かい話はあとでちゃんと書きます。）最終的にあるスレーター行列式への、1体演算子の総和で書かれた多体ハミルトニアンの作用は、

$$
\begin{align*}
&\mathcal{H}\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|\\

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

となります。ここで前半の

$$
\sum_i \hat{H}_i
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

と、最後の

$$
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

$$

を見比べてみると、ハミルトニアンを作用させることで、「スレーター行列式の中の$\varphi_\lambda$を$\varphi_\delta$に変えたもの（の$\delta$の総和）と、$\varphi_\mu$を$\varphi_\delta$に変えたもの（の$\delta$の総和）と$\cdots$という総和になっていることがわかります。そこで、
- スレーター行列式の中の状態$\varphi_\kappa$を「消す」演算子$\hat{a}_\kappa$
- スレーター行列式に状態$\varphi_\delta$を「加える」演算子$\hat{a}_\delta^\dagger$
  
を上手く考えると上記の作用が

$$
\sum_i \hat{H}_i
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|

=
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

のような感じで書けんじゃん！というのが、ざっくりとしたハミルトニアンの第二量子化表示の概要です。なんかいけそうですね。


導入はここまでとして、上記のようなゴールを念頭に、ここからもうちょい真面目に書いて行きます。

余談ですが、私は正直なところ、第二量子化を雰囲気で使っていたので、「電子を生成する・付け加える演算子？消す演算子？電子ってそんな空間にポンポン付け加えたり消したりできるもんなの？」みたいな疑問をずっと持ちながら生き長らえてていたのですが、改めて第二量子化をちゃんと勉強したら「電子を生成/消滅するというか、スレーター行列式に加える・消す演算子なのか」と納得することができました。^[この文を公開するのはかなり恥ずかしいのですが、多分私と同じようなことで苦悶している方も数人程度はいるだろうと、恥を忍んで（脚注にも逃げず）本文に書いておきます。]^[とはいえ、結局スレーター行列式に/から波動関数を付け加える/消す操作は「空間に電子を加える/消す」ことを意味するのですが、当面は空間内の電子数は一定の状況を考えるので、基本的には電子を消す「消滅演算子」$\hat{a}_\delta$と電子を付け加える「生成演算子」$\hat{a}^\dagger_\kappa$はセットで現れ、二つ合わせて「電子の状態を$\delta$から$\kappa$に変える演算子$\hat{a}^\dagger_\kappa\hat{a}_\delta$のように使われることになります。もしかしたらさらに進んでいくと本当に「付け加えたり」「消したり」することもあるのでしょうが（例えば電子ではないですが、格子振動を生成消滅演算子で表した場合は、あるエネルギーの振動（フォノン）を発生させることが、生成演算子を作用させることに対応してたりしたような気がします）]



## 多体電子系のシュレーディンガー方程式

最終的な目標が見えたところで、改めて、上記のような表現を実現するためにどのような考え方をしているのか、そして演算子（生成演算子・消滅演算子）の定義や、性質を整理していきます。最後に、そのように定義することで実際に先ほど期待したような第二量子化の表現が実現できることを見ていきます。

初めに、多体電子系のシュレーディンガー方程式から出発します。シュレーディンガー方程式は多体の波動関数$\Phi(\tau_1,\tau_2,\cdots,\tau_N)$とハミルトニアン$\mathcal{H}$を用いて、

$$
\mathcal{H}\Phi(\tau_1,\tau_2,\cdots,\tau_N)
=
E
\Phi(\tau_1,\tau_2,\cdots,\tau_N)
$$

と書かれます。ここで$\tau_i$は、$i$番目の電子の位置座標$\boldsymbol{r}$とスピン座標$\sigma=\uparrow, \downarrow$を合わせた変数です。

ここで多体の波動関数$\Phi(\tau_1,\tau_2,\cdots,\tau_N)$は一般には1変数関数の積に分解できず、多変数関数となってしまいます。

しかしここで[この章](https://zenn.dev/ponzumai/articles/tight-binding-model-spin)で導入した反対称化された波動関数であるスレーター行列式を、完全正規直交性を満たす軌道部分の波動関数$\varphi_a, \varphi_b\cdots$と、上向き/下向きスピン状態を合わせたもの；$\varphi_{a\uparrow}(\boldsymbol{r}), \varphi_{a\downarrow}(\boldsymbol{r}),\varphi_{b\uparrow}(\boldsymbol{r}), \varphi_{b\downarrow}(\boldsymbol{r})\cdots$、または別の表現をするとスピン関数$\alpha(\sigma)$、$\beta(\sigma)$との積$\varphi_a(\boldsymbol{r})\alpha(\sigma), \varphi_a(\boldsymbol{r})\beta(\sigma),\cdots$を考え、軌道部分の関数をラベルする量子数（例えば水素様原子中の電子の場合だと$(n,l,m)$）とスピン状態のラベル$\uparrow,\downarrow$をまとめて、ギリシャ文字$\lambda,\mu,\xi$などと表した「スピン軌道関数」$\varphi_\lambda(\tau_1), \varphi_\mu(\tau_2),\cdots,\varphi_\xi(\tau_N)$を用いて以下のように：

$$
\begin{align*}
 \left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right| &\equiv
\frac{1}{\sqrt{N!}}
\begin{vmatrix}
\varphi_\lambda(\tau_1) & \varphi_\mu(\tau_1) & \cdots & \varphi_\xi(\tau_1)\\
\varphi_\lambda(\tau_2) & \varphi_\mu(\tau_2) & \cdots & \varphi_\xi(\tau_2)\\
& \cdots & \cdots\\
\varphi_\lambda(\tau_N) & \varphi_\mu(\tau_N) & \cdots & \varphi_\xi(\tau_N)
\end{vmatrix}\\

&=\frac{1}{\sqrt{N!}}\sum_P(-1)^P\hat{P}\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N})\\

\end{align*}
$$

と定義すると、（ここで右辺二つ目の式は、行列式を、座標を入れ替える置換演算子$\hat{P}$と、置換の回数だけ$(-1)$を掛ける$(-1)^P$で表したものです。）任意の多体の波動関数$\Phi(\tau_1,\tau_2,\cdots,\tau_N)$をスレーター行列式の線形結合として

$$
\Phi(\tau_1,\tau_2,\cdots,\tau_N)
=
\sum_{\lambda,\mu,\cdots,\xi}C_{\lambda,\mu,\cdots\xi}
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

と展開することができます。

これは[格子ベクトルと逆格子ベクトル](https://zenn.dev/ponzumai/articles/tight-binding-model-lattice-vec)の章で多変数関数を1変数のFourier展開の積で書き直したようなことを、さらに多くの変数とFourier級数に限らない完全正規直交関数系でやってるようなもので、新しいのはフェルミオン（電子）の波動関数が満たすべき反対称性を満たすためにスレーター行列式の和になっているくらいで、そこまでややこしいことはしていません。

さて、多体の波動関数を上記のように完全正規直交系を用いて展開したという前提の元、話を進めていきましょう。このように波動関数を展開する関数系を指定したとき、その関数系を**基底**と呼んだりします。特に完全正規直交系であるとき、**完全正規直交基底**と呼んだりします。

このように基底を指定したとき、シュレーディンガー方程式の中のハミルトニアン演算子$\mathcal{H}$は、「任意の完全正規直交系$\{\varphi_\nu\}$で書かれたスレーター行列式に作用した結果どうなるか？」という演算子を考えると書き換えることができそうです。

「は？」となりそうなので、一つ具体例を挟んでおきます。1電子ハミルトニアン$\hat{H}=-\hbar^2/2m \nabla^2 + U(\boldsymbol{r})$でのシュレーディンガー方程式はお馴染み$\frac{-\hbar^2}{2m}\nabla^2\varphi(\boldsymbol{r}) = \varepsilon\varphi(\boldsymbol{r})$ですが、ここで波動関数を平面波（これは波動関数の境界条件に周期的境界条件を採用すると完全正規直交性を満たす関数系となっています）$e^{i\boldsymbol{k}\cdot\boldsymbol{r}}$で$\varphi(\boldsymbol{r}) = \sum_{\boldsymbol{k}}C_{\boldsymbol{k}}e^{i\boldsymbol{k}\cdot\boldsymbol{r}}$と展開した場合を考えると、任意の平面波$e^{i\boldsymbol{k}\cdot\boldsymbol{r}}$に対して$\hat{H}e^{i\boldsymbol{k}\cdot\boldsymbol{r}} = \{(\hbar^2/2m)k^2 + U(\boldsymbol{r})\}e^{i\boldsymbol{k}\cdot\boldsymbol{r}}$となるので、波動関数を平面波で展開している、という状況が続く限りは、「平面波の波数の絶対値を抜き出す演算子」$\hat{k}$（注：これは私が勝手に定義した演算子で、こんなものは普通ありません）を定義してハミルトニアンを$\hat{H} = (\hbar^2/2m )\hat{k}^2 + U(\boldsymbol{r})$と書いても問題ないわけです。なんか余計混乱しそうで蛇足だったかもしれませんがそういうことです。

さて、話を戻して、多体の波動関数を

$$
\Phi(\tau_1,\tau_2,\cdots,\tau_N)=\sum_{\lambda,\mu,\cdots,\xi}C_{\lambda,\mu,\cdots\xi}
\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|
$$

と展開したとき、その各展開項の、完全正規直交性を満たす関数系$\{\varphi_\nu\}$で構成されるスレーター行列式$\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|$に対して多体のハミルトニアン$\mathcal{H}$を作用させたとき一般にどのようなことがおきるか？

$$
\mathcal{H}\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right| = ??

$$

を考えることで、多体のハミルトニアンを何かしらの演算子を定義して書き換えることができそうです。これが結局、第二量子化表示と呼ばれるものです。




## スレーター行列式への1体演算子の作用

「何かしらの演算子」とは何かというと、それがまさしく冒頭で述べた「生成・消滅演算子」と呼ばれるものです。正規直交関数系$\{\varphi_\nu\}$を指定したときに、任意のスレーター行列式$\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|$に対して多体のハミルトニアンを作用させたとして、何かしら一般的に「こうなる」と言えるものはあるのでしょうか。

ここからは冒頭と重なる部分もありますが、改めて1から書いていきます。

まず、多体のハミルトニアン$\mathcal{H}$が1電子のハミルトニアン$\hat{H}$の和で

$$
\mathcal{H} = \sum_i\hat{H}_i
$$

と書かれている場合について考えます（というか本稿ではその場合のみ考えて、2体演算子の場合は扱いません）。

このとき多体のハミルトニアンをスレーター行列式に作用させてみると、冒頭で述べた通り

$$
\begin{align*}
&\mathcal{H}\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|\\



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
    \right|

\end{align*}
$$

と、スレーター行列式の各1電子の波動関数のどれか一つに対して、1体ハミルトニアンが作用したスレーター行列式の総和となります。

証明は以下のようにして、実際にスレーター行列式に作用させた後に、出てきた項を並べ替えることで示せます。（証明というか確認みたいな感じですが）

$$
\begin{align*}
&\mathcal{H}\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|\\

&=
\sum_i\hat{H}_i
\frac{1}{\sqrt{N!}}\sum_P(-1)^P\hat{P}\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N})\\

&=
\frac{1}{\sqrt{N!}}
\left[
\left\{
    \hat{H}_1+ \hat{H}_2 + \cdots\hat{H}_N
    \right\}
    \left\{
(-1)^0\varphi_\lambda(\tau_1)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N})
        \right. \right.\\
        
&\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>

+(-1)^1\varphi_\lambda(\tau_2)\varphi_\mu(\tau_1)\cdots\varphi_\xi(\tau_{N})
+\cdots
        \left.\left.
        
        \right\}
    \right]\\

&=
\frac{1}{\sqrt{N!}}
\left[
    (-1)^0\left(\hat{H}_1\varphi_\lambda(\tau_1)\right)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N}) 
    +
    (-1)^0\varphi_\lambda(\tau_1)\left(\hat{H}_2\varphi_\mu(\tau_2)\right)\cdots\varphi_\xi(\tau_{N}) 
    +\cdots
    \right.
    \\
&\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>\>
\left.
+
(-1)^1\varphi_\lambda(\tau_2)\left(\hat{H}_1\varphi_\mu(\tau_1)\right)\cdots\varphi_\xi(\tau_{N})
+
(-1)^1\left(\hat{H}_2\varphi_\lambda(\tau_2)\right)\varphi_\mu(\tau_1)\cdots\varphi_\xi(\tau_{N})
+\cdots
\right]\\

&=
\frac{1}{\sqrt{N!}}
\left[
    (-1)^0\left(\hat{H}_1\varphi_\lambda(\tau_1)\right)\varphi_\mu(\tau_2)\cdots\varphi_\xi(\tau_{N}) 
    +
    (-1)^1\left(\hat{H}_2\varphi_\lambda(\tau_2)\right)\varphi_\mu(\tau_1)\cdots\varphi_\xi(\tau_{N})
    + \cdots\right]
\\

&\>\>\>\>+
\frac{1}{\sqrt{N!}}
\left[

    (-1)^0\varphi_\lambda(\tau_1)\left(\hat{H}_2\varphi_\mu(\tau_2)\right)\cdots\varphi_\xi(\tau_{N}) 
    
+
(-1)^1\varphi_\lambda(\tau_2)\left(\hat{H}_1\varphi_\mu(\tau_1)\right)\cdots\varphi_\xi(\tau_{N})
+

+\cdots
\right]\\

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
    \right|.

\end{align*}
$$


ここで、繰り返しになりますが、1電子波動関数$\varphi_\lambda, \varphi_\mu,\cdots\varphi_\xi$は完全系であると仮定したので、それぞれに1体のハミルトニアンが作用した後の関数を$\hat{H}\varphi_\kappa(\tau) = \sum_{\delta}c_{\delta}^{\kappa}\varphi_{\delta}(\tau)$と展開でき、また正規直交性より、展開係数

$$
c_{\delta}^{\kappa} = \int \varphi_{\delta}^*(\tau)\hat{H}\varphi_{\kappa}(\tau)d\tau
=
\sum_{\sigma=\uparrow,\downarrow}\int \varphi_{d}^*(\boldsymbol{r})\gamma^*_d(\sigma)\hat{H}\varphi_{k}(\boldsymbol{r})\gamma_k(\sigma)d\boldsymbol{r}
$$

となります。ここで最後の式は$\varphi_\delta(\tau) = \varphi_d(\boldsymbol{r})\gamma_d(\sigma), \gamma(\sigma) = \alpha(\sigma), \beta(\sigma)$として$\tau$の積分を書き換えたものです。

これを$\int \varphi_{\delta}^*(\tau)\hat{H}\varphi_{\kappa}(\tau)d\tau\equiv H_{\delta\kappa}$と書くと、最終的にあるスレーター行列式への、1体演算子の総和で書かれた多体ハミルトニアンの作用は、

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

となります。これはハミルトニアンを作用させることで、「スレーター行列式の中の$\varphi_\lambda$を$\varphi_\delta$に変えたもの（の$\delta$の総和）と、$\varphi_\mu$を$\varphi_\delta$に変えたもの（の$\delta$の総和）と$\cdots$という総和になっていることがわかります。そこで、繰り返しになりますが、
- スレーター行列式の中の状態$\varphi_\kappa$を「消す」演算子$\hat{a}_\kappa$
- スレーター行列式に状態$\varphi_\delta$を「加える」演算子$\hat{a}_\delta^\dagger$
  
を上手く考えることで上記の作用を

$$
\sum_i \hat{H}_i
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|

=
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

のように書き換えたいわけです。


## 生成・消滅演算子の導入とその満たすべき性質

さて、
- スレーター行列式の中の状態$\varphi_\kappa$を「消す」演算子$\hat{a}_\kappa$
- スレーター行列式に状態$\varphi_\delta$を「加える」演算子$\hat{a}_\delta^\dagger$

はどのように定義すると上手くいくのかをここから見ていきます。

### 生成演算子の定義

まず、スレーター行列式にラベル$\delta$（これは軌道関数のラベルと、スピン関数の両方を合わせたラベルとします）で示される1電子状態$\varphi_\delta$を付け加える**生成演算子**を$\hat{a}_\delta^\dagger$と定義します。特にスレーター行列式において行列式の性質から、

$$
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
=
-
\left|
\varphi_\mu\varphi_\lambda\cdots\varphi_\xi
    \right|
$$

のように、1電子状態の順序を変えると符号が変わってしまいますので、どの場所に付け加えるか決めておく必要があります。そこで生成演算子$\hat{a}_\delta^\dagger$は、

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

と、スレーター行列式の**一番手前**に状態$\varphi_\delta$を付け加えるものとします。^[上記のように生成演算子を定義すると、$\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|$のようなスレーター行列式を、生成演算子を用いて書くことができそうです。特に、電子を一つも含まない「状態」を$\ket{0}$と書いて、$\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\xi\right|=\hat{a}_\lambda^\dagger\hat{a}_\mu^\dagger\cdots\hat{a}_\xi^\dagger\ket{0}$と書くことができそうです。また、あるスレーター行列式を、「$\varphi_\lambda$状態の電子が1つ、$\varphi_\mu$状態の電子が1つ$\cdots$」のようにある状態を取る電子の数で表すことができそうです。これを数表示と呼んだりします。$\cdots$が、上式からわかるようにこのような話を進めていくには波動関数のベクトル表現、そしてブラケット記法が必要になり、さらに多体電子系のベクトル空間（Hilbert空間）は1電子系のHilbert空間のテンソル積であり…等とやらないといけないことがたくさん出てくるので、本稿では扱わず、電子を一つも含まない状態も雑な表現でごまかしています。ちゃんと勉強したい方は量子力学の教科書を見てください。（いつか心身に余裕ができたらこのあたりもちゃんと書いていきたいところですが）]

また、この時状態$\varphi_\delta$を付け加えると同時に、スレーター行列式の規格化係数$1/\sqrt{N!}$を$1/\sqrt{(N+1)!}$に変えているものとします。（今回の書き方では規格化係数は明示的に書かれないので、ぱっと見はわからないのですが）


また、スレーター行列式の定義から、既にスレーター行列式の中のどこかに状態$\varphi_\delta$があった場合、生成演算子$\hat{a}_\delta^\dagger$を作用させた結果は$0$となります：

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

続いて先ほどハミルトニアンを書き換えたように、スレーター行列式から状態を消す**消滅演算子**を定義する必要もあります。ひとまず、記号としては状態$\varphi_\kappa$をスレーター行列式から「消す」演算子を$\hat{a}_\kappa$と定義します。ただ、こちらも

$$
\hat{a}_\kappa|\cdots \varphi_\kappa\cdots|=|\cdots \cdots|(?)
$$

のように好き勝手に「消す」と、並べ替えに対する符号の変化と整合性が取れなくなってしまいます。例えば、生成演算子を二つ作用させたスレーター行列式

$$
\hat{a}_\alpha^\dagger\hat{a}_\beta^\dagger|\cdots| = |\varphi_\alpha\varphi_\beta\cdots|
$$

に、「好きに消せる」消滅演算子（？）$\hat{a}_\alpha, \hat{a}_\beta$を作用させてみます

$$
\hat{a}_\alpha\hat{a}_\beta |\varphi_\alpha\varphi_\beta\cdots| = |\cdots|.
$$

ところで、スレーター行列式は順番を入れ替えると符号が変わるので、$|\varphi_\alpha\varphi_\beta\cdots| = -|\varphi_\beta\varphi_\alpha\cdots|$です。こちらにも同じく消滅演算子（？）を作用させてみると、

$$
\hat{a}_\alpha\hat{a}_\beta (-1)|\varphi_\beta\varphi_\alpha\cdots| = - |\cdots|.
$$

以上より、

$$
\hat{a}_\alpha\hat{a}_\beta |\varphi_\alpha\varphi_\beta\cdots| =
\hat{a}_\alpha\hat{a}_\beta (-1)|\varphi_\beta\varphi_\alpha\cdots| \\
\Rightarrow
|\cdots| = - |\cdots|(?)
$$

と、変な結論になってしまいます。というわけで消滅演算子も、状態を「消す」時に何かしらの順番を考慮する必要がありそうです。（ここでも、生成演算子の場合と同じく規格化係数を$1/\sqrt{N!}$から$1/\sqrt{(N-1)!}$に変えているものとします）

これは状態$\kappa$についての消滅演算子$\hat{a}_\kappa$を、「スレーター行列式の先頭にある1電子状態$\varphi_\kappa$を消す、先頭になければ先頭まで交換して移動してから消す」と定義すると上手くいきます。つまり、

$$
\hat{a}_\kappa|\varphi_\kappa\cdots \cdots|=|\cdots \cdots|
$$

または、

$$
\hat{a}_\kappa|(n個の状態)\varphi_\kappa \cdots|=
\hat{a}_\kappa(-1)^n|\varphi_\kappa(n個の状態) \cdots|=(-1)^n|\cdots \cdots|
$$

と定義します。ややまどろっこしい定義に思えますが、これは生成演算子との関係

$$
\hat{a}_\kappa\hat{a}_\kappa^\dagger = 1
$$

として定義しているのだと考えるとやや腑に落ちるかと思います。（次に見る「反交換関係」を整理していくことで、これが上手い定義になっているのだとわかります。）

また、スレーター行列式に消滅演算子$\hat{a}_\kappa$の示す状態$\varphi_\kappa$が含まれない場合は、消滅演算子を作用させた結果は$0$を返すものと定義します。：

$$
\hat{a}_\kappa|\cdots(\varphi_\kappa 無し) \cdots|=0.
$$

### 生成消滅演算子の満たす性質：反交換関係

上記のように生成・消滅演算子を定義すると、スレーター行列式の満たすべき性質を利用して、以下のような生成・消滅演算子間の交換の規則「反交換関係」を導出することができます。^[これはスレーター行列式を上手く再現するための「定義」のように見えなくもないのですが、生成消滅演算子の定義部分でスレーター行列式と関連付けたので、やっぱりそこから導出されるものであり定義ではないのでしょう。まああまり細かいことは気にせず行きましょう]

:::message
生成消滅演算子の満たす反交換関係

$$
\hat{a}_\kappa\hat{a}_\zeta + \hat{a}_\zeta\hat{a}_\kappa  
= 
\hat{a}_\delta^\dagger\hat{a}_\eta^\dagger + \hat{a}_\eta^\dagger\hat{a}_\delta^\dagger=0\\
\hat{a}_\delta^\dagger\hat{a}_\kappa + \hat{a}_\kappa\hat{a}_\delta^\dagger=
\delta_{\delta\kappa}
$$

:::

#### 消滅演算子の反交換関係

まずは順番に、消滅演算子の反交換関係

$$
\hat{a}_\kappa\hat{a}_\zeta + \hat{a}_\zeta\hat{a}_\kappa  =0
$$

を示していきます。消滅演算子以外の反交換関係も含め、スレーター行列式への作用を考えることで示せます。

まず作用する先のスレーター行列式に、二つの消滅演算子に対応する状態がない場合は、演算子を作用させた結果が$0$になり、等式は成り立ちます。

$$
(\hat{a}_\kappa\hat{a}_\zeta + \hat{a}_\zeta\hat{a}_\kappa )
|\varphi_\lambda\cdots (\varphi_\kappa,\varphi_\zeta どちらか、又はどちらも無い)\cdots|

=0.
$$

そこでスレーター行列式に$\varphi_\kappa,\varphi_\zeta$どちらも含まれる場合を考え、$\varphi_\kappa$が$\varphi_\zeta$よりも左にある場合と、右にある場合で場合分けします。

$\varphi_\kappa$が$\varphi_\zeta$よりも左にある場合

$$
\begin{align*}
&\hat{a}_\kappa\hat{a}_\zeta
|\cdots(n_1 個) \cdots \varphi_\kappa\cdots (n_2個) \cdots \varphi_\zeta\cdots|
\\

&=
\hat{a}_\kappa\hat{a}_\zeta
(-1)^{n_1}(-1)^{n_1+n_2+1}|\varphi_\kappa\varphi_\zeta\cdots(n_1 個) \cdots \cdots (n_2個) \cdots \cdots|\\

&=
(-1)^{2n_1+n_2+1}|\cdots \cdots  \cdots|,\\

&\hat{a}_\zeta\hat{a}_\kappa
|\cdots(n_1 個) \cdots \varphi_\kappa\cdots (n_2個) \cdots \varphi_\zeta\cdots|
\\

&=
\hat{a}_\zeta\hat{a}_\kappa
(-1)^{n_1}(-1)^{n_1+n_2}|\varphi_\zeta\varphi_\kappa\cdots(n_1 個) \cdots \cdots (n_2個) \cdots \cdots|\\

&=
(-1)^{2n_1+n_2}|\cdots\cdots \cdots|.

\end{align*}
$$

一方$\varphi_\kappa$が$\varphi_\zeta$よりも右にある場合も、

$$
\begin{align*}
&\hat{a}_\kappa\hat{a}_\zeta
|\cdots(n_1 個) \cdots \varphi_\zeta\cdots (n_2個) \cdots \varphi_\kappa\cdots|
\\

&=
\hat{a}_\kappa\hat{a}_\zeta
(-1)^{n_1}(-1)^{n_1+n_2}|\varphi_\zeta\varphi_\kappa\cdots(n_1 個) \cdots \cdots (n_2個) \cdots \cdots|\\

&=
(-1)^{2n_1+n_2}|\cdots \cdots  \cdots|,\\

&\hat{a}_\zeta\hat{a}_\kappa
|\cdots(n_1 個) \cdots \varphi_\zeta\cdots (n_2個) \cdots \varphi_\kappa\cdots|
\\

&=
\hat{a}_\zeta\hat{a}_\kappa
(-1)^{n_1}(-1)^{n_1+n_2 + 1}|\varphi_\kappa\varphi_\zeta\cdots(n_1 個) \cdots \cdots (n_2個) \cdots \cdots|\\

&=
(-1)^{2n_1+n_2 + 1}|\cdots\cdots \cdots|.

\end{align*}
$$


のように、$\varphi_\kappa$、$\varphi_\zeta$の順序によらず、消滅演算子の順序を入れ替えるとその結果の符号が逆になるため、任意のスレーター行列式に対して

$$
(\hat{a}_\kappa\hat{a}_\zeta + \hat{a}_\zeta\hat{a}_\kappa)|\cdots| = 0
$$

となることが示せました。

#### 生成演算子の反交換関係

続いて生成演算子の方の関係ですが、こちらは作用先のスレーター行列式に、二つある演算子のラベルのどちらか（に対応する1電子状態）があれば自動的に$0$となります。一方どちらも含まない場合を考えると、

$$
\hat{a}_\delta^\dagger\hat{a}_\eta^\dagger |\cdots| = 
|\varphi_\delta\varphi_\eta\cdots|,\\

\hat{a}_\eta^\dagger\hat{a}_\delta^\dagger|\cdots| = 
|\varphi_\eta\varphi_\delta\cdots|
=
-|\varphi_\delta\varphi_\eta\cdots|
$$

より、こちらは比較的シンプルに、任意のスレーター行列式に対して$(\hat{a}_\delta^\dagger\hat{a}_\eta^\dagger + \hat{a}_\eta^\dagger\hat{a}_\delta^\dagger)|\cdots|=0$を示せました。

#### 生成演算子と消滅演算子の交換関係

最後に生成演算子と消滅演算子の間の交換関係

$$
\hat{a}_\delta^\dagger\hat{a}_\kappa + \hat{a}_\kappa\hat{a}_\delta^\dagger=
\delta_{\delta\kappa}
$$


について、二つの演算子のラベルが異なる場合と同じ場合で場合分けして考えましょう。

ラベルが異なる場合

演算子$\hat{a}_\delta^\dagger, \hat{a}_\kappa$を考えるとして、これまでと似たように、作用させる先のスレーター行列式に$\varphi_\delta$がなく、かつ$\varphi_\kappa$がある場合以外はゼロになるので、ゼロにならない場合のみ考えます。

$$
\hat{a}_\delta^\dagger\hat{a}_\kappa
|\cdots (n個) \varphi_\kappa\cdots|
=
\hat{a}_\delta^\dagger \hat{a}_\kappa
(-1)^n
| \varphi_\kappa\cdots (n個)\cdots|
=
(-1)^n
| \varphi_\delta\cdots (n個)\cdots|,\\

\hat{a}_\kappa\hat{a}_\delta^\dagger
|\cdots (n個) \varphi_\kappa\cdots|
=
\hat{a}_\kappa
| \varphi_\delta\cdots (n個)\varphi_\kappa\cdots|
=
(-1)^{n+1}
| \varphi_\delta\cdots (n個)\cdots|
$$

と、符号が逆になるので、ラベルが異なる場合、交換関係$\hat{a}_\delta^\dagger\hat{a}_\kappa + \hat{a}_\kappa\hat{a}_\delta^\dagger = 0$を示せます。

一方ラベルが同じ場合（$\delta$とします）は、ちょっとややこしいのですが、まずスレーター行列式に$\varphi_\delta$がある場合：

「消滅」「生成」の順番

$$
\hat{a}_\delta\hat{a}_\delta^\dagger
|\cdots (n個) \varphi_\delta\cdots|
=0
$$

はゼロになります。一方、「生成」「消滅」の順番は、

$$
\begin{align*}
\hat{a}_\delta^\dagger\hat{a}_\delta
|\cdots (n個) \varphi_\delta\cdots|
&=
\hat{a}_\delta^\dagger \hat{a}_\delta
(-1)^n
| \varphi_\delta\cdots (n個)\cdots|\\
&=
(-1)^n
| \varphi_\delta\cdots (n個)\cdots|\\

&=
(-1)^n(-1)^n
| \cdots (n個)\varphi_\delta\cdots|\\

&=
|\cdots (n個) \varphi_\delta\cdots|
\end{align*}
$$

と、スレーター行列式は演算子の積$\hat{a}_\delta^\dagger\hat{a}_\delta$の固有状態になっていることがわかり、その場合$\hat{a}_\delta^\dagger\hat{a}_\delta + \hat{a}_\delta\hat{a}_\delta^\dagger  = 1$を満たすことがわかります。

一方まずスレーター行列式に$\varphi_\delta$が無い場合：

今度は「生成」「消滅」の順番は、

$$
\hat{a}_\delta^\dagger\hat{a}_\delta
|\cdots ( \varphi_\delta 無し)\cdots|
=0
$$

ゼロになります。一方「消滅」「生成」の順番は、

$$
\begin{align*}
\hat{a}_\delta\hat{a}_\delta^\dagger
|\cdots ( \varphi_\delta 無し)\cdots|
&=
\hat{a}_\delta
| \varphi_\delta\cdots\cdots|\\
&=

| \cdots ( \varphi_\delta 無し)\cdots|\\

\end{align*}
$$

と、今度はこちらが固有状態になっており、ここでも$\hat{a}_\delta^\dagger\hat{a}_\delta + \hat{a}_\delta\hat{a}_\delta^\dagger  = 1$$を満たすことがわかります。

以上から、任意のスレーター行列式に対して

$$
\hat{a}_\delta^\dagger\hat{a}_\kappa + \hat{a}_\kappa\hat{a}_\delta^\dagger=
\delta_{\delta\kappa}
$$

が成り立つことが示せました。なおこれは、消滅演算子の状態の「消し方」$\hat{a}_\delta\hat{a}_\delta^\dagger = 1$の一般的な場合とも見れます。


## 1体演算子のみの多体ハミルトニアンを第二量子化表示で表す

さて、以上見てきた生成・消滅演算子を用いて、冒頭で書いたような一体演算子の総和で書かれた多体ハミルトニアンのスレーター行列式への作用を書き換えるのが（ざっくりとした）第二量子化の考え方です。（本当は場の演算子を使って考えることを「第二量子化」と呼ぶっぽいような気もするのですが、生成・消滅演算子で多体ハミルトニアンを書き換えることを「第二量子化」と呼んでいる事例も多いので、そのあたりは適当にやっておきます）

1体演算子の総和の、スレーター行列式への作用を再掲しておくと、

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

でした。ここまでの生成・消滅演算子の定義や性質を踏まえると、上式を生成・消滅演算子を用いて書き換えることは造作ないように思えますね。早速やっていきましょう。


まず3式目と4式目から、適当な項をとりだしてみます：

$$
\left|
\varphi_\lambda\varphi_\mu\cdots\left(\hat{H}\varphi_\kappa\right)\cdots\varphi_\xi
    \right|

=

\sum_{\delta}H_{\delta\kappa}\left|
    \varphi_{\lambda}\varphi_\mu\cdots\varphi_\delta\cdots\varphi_\xi
    \right|.
$$

これは、「ハミルトニアンを作用させた結果、もとのスレーター行列式$\left|\varphi_\lambda\varphi_\mu\cdots\varphi_\kappa\cdots\varphi_\xi \right|$で$\varphi_\kappa$があった場所の状態を$\varphi_\delta$に変えて、そののち$\int \varphi_\delta^*(\tau)\hat{H}\varphi_\kappa(\tau)d\tau\equiv H_{\delta\kappa}$を掛けて、かつ全ての$\delta$について足し合わせる」ことを意味しています。

あるスレーター行列式の$\varphi_\kappa$を、場所を変えずに$\varphi_\delta$に変えるには、以下のように演算子$\hat{a}^\dagger_\delta\hat{a}_\kappa$をかければ良いことがわかります：

$$
\begin{align*}
\hat{a}^\dagger_\delta\hat{a}_\kappa
\left|\varphi_\lambda\varphi_\mu\cdots(n個)\varphi_\kappa\cdots\varphi_\xi \right|
&=
\hat{a}^\dagger_\delta\hat{a}_\kappa
(-1)^n
\left|\varphi_\kappa\varphi_\lambda\varphi_\mu\cdots(n個)\cdots\varphi_\xi \right|
\\

&=

(-1)^n
\left|\varphi_\delta\varphi_\lambda\varphi_\mu\cdots(n個)\cdots\varphi_\xi \right|
\\

&=

(-1)^n(-1)^n
\left|\varphi_\lambda\varphi_\mu\cdots(n個)\varphi_\delta\cdots\varphi_\xi \right|
\\

&=
\left|\varphi_\lambda\varphi_\mu\cdots(n個)\varphi_\delta\cdots\varphi_\xi \right|.
\end{align*}
$$

これは「$\delta$に変える」対象の状態が何であろうと、符号などを心配せずに共通して$\hat{a}_\delta^\dagger\hat{a}_\kappa$をかければOKであることを示しています。

従って、当初の期待通り、「$\delta$に変える元の全ての状態」についての和と、「変えた後の$\delta$について全ての状態」についての和を取って、

$$
\begin{align*}
\mathcal{H}\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
&=\sum_i \hat{H}_i
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|\\

&=
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|

\end{align*}
$$

$$
H_{\delta\kappa}\equiv \int \varphi_{\delta}^*(\tau)\hat{H}\varphi_{\kappa}(\tau)d\tau

$$

と書けることが示せました！

なお、ここで当初は予想していなかった生成（消滅）演算子の利点として、「生成（消滅）演算子の作用先のスレーター行列式にラベルと同じ状態があれば（無ければ）ゼロを返す」という性質から、スレーター行列式の粒子数（電子数）や、元々含まれていた状態を気にせず、「消す状態のラベル$\kappa$」は全ての状態について和を取ればよいことが分かります。

ここで一旦まとめると、ハミルトニアンが1体演算子の総和で書かれた多体のシュレーディンガー方程式


$$
\mathcal{H}\Phi(\tau_1,\tau_2,\cdots,\tau_N)
=
E
\Phi(\tau_1,\tau_2,\cdots,\tau_N),\\

\mathcal{H} = \sum_i\hat{H}_i
$$

において、完全正規直交系$\{\varphi_\nu\}$を用いて多体の波動関数$\Phi(\tau_1,\tau_2,\cdots,\tau_N)$をスレーター行列式の線形結合として


$$
\Phi(\tau_1,\tau_2,\cdots,\tau_N)
=
\sum_{\lambda,\mu,\cdots,\xi}C_{\lambda,\mu,\cdots\xi}
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
$$

と展開したとき、このような完全正規直交系$\{\varphi_\nu\}$を選んでいるという前提のもと、ハミルトニアンの任意のスレーター行列式への作用

$$
\begin{align*}
\mathcal{H}\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|
&=\sum_i \hat{H}_i
\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|\\

&=
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa\left|
\varphi_\lambda\varphi_\mu\cdots\varphi_\xi
    \right|

\end{align*}
$$

$$
H_{\delta\kappa}\equiv \int \varphi_{\delta}^*(\tau)\hat{H}\varphi_{\kappa}(\tau)d\tau

$$

から、1体演算子の総和で書かれた多体のハミルトニアンは

$$
\mathcal{H} = \sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa
$$

と書けることが分かりました。

## ハミルトニアンがスピンについての項を含まない場合

なお、1体のハミルトニアンが例えば$\hat{H} = -(\hbar^2/2m)\nabla^2$のようにスピンについての項を含まない場合、 $\varphi_{\delta}(\tau) = \varphi_d(\boldsymbol{r})\gamma_d(\sigma), \varphi_{\kappa}(\tau) = \varphi_k(\boldsymbol{r})\gamma_k(\sigma)$、$\gamma_d(\sigma), \gamma_k(\sigma) = \alpha(\sigma)$または$\beta(\sigma)$のどちらか、として、

$$
\begin{align*}
H_{\delta\kappa}&\equiv \int \varphi_{\delta}^*(\tau)\hat{H}\varphi_{\kappa}(\tau)d\tau\\
&=
\sum_\sigma\int \varphi_{d}^*(\boldsymbol{r})\gamma_d(\sigma)\hat{H}\varphi_{k}(\boldsymbol{r})\gamma_k(\sigma)d\boldsymbol{r}\\

&=
\sum_\sigma\gamma_d(\sigma)\gamma_k(\sigma)\int \varphi_{d}^*(\boldsymbol{r})\hat{H}\varphi_{k}(\boldsymbol{r})d\boldsymbol{r}\\

&=
\delta_{\gamma_d,\gamma_k}
\int \varphi_{d}^*(\boldsymbol{r})\hat{H}\varphi_{k}(\boldsymbol{r})d\boldsymbol{r}
\end{align*}
$$

となり、

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

と、スピン部分については同じスピン状態を取り、軌道関数でハミルトニアンを挟んだ積分を係数とする第二量子化表示でハミルトニアンを書くことができます。

ただし、2つめの式で出てきた$\hat{a}^\dagger_{d,\gamma}$は、スピン軌道関数のうち軌道関数部分が$\varphi_d(\boldsymbol{r})$、スピン関数部分がそれぞれ$\gamma=\alpha$なら$\alpha(\sigma)$となるスピン軌道関数$\varphi_\delta(\tau)=\varphi_d(\boldsymbol{r})\alpha(\sigma)$、$\gamma=\beta$なら$\beta(\sigma)$となるスピン軌道関数$\varphi_\delta(\tau)=\varphi_d(\boldsymbol{r})\beta(\sigma)$を生成する演算子です。（消滅演算子$\hat{a}_{k,\gamma}$も同様。）

なお、スピン部分のラベルは、$\sigma$を用いて$\sum_\sigma\hat{a}_{d,\sigma}^\dagger\hat{a}_{k,\sigma}$のように書かれることも多いです。

また、$H_{dk}$は軌道関数のみで挟んだハミルトニアンの積分$H_{dk}\equiv\int \varphi_{d}^*(\boldsymbol{r})\hat{H}\varphi_{k}(\boldsymbol{r})d\boldsymbol{r}$として定義しました。

## 多体のシュレーディンガー方程式を解く

さて、以上のように、1体演算子の総和で書かれた多体のハミルトニアンを

$$
\begin{align*}
\mathcal{H}
&=
\sum_{\delta\kappa}\int\varphi_\delta(\tau)\hat{H}\varphi_\kappa(\tau)d\tau\hat{a}^\dagger_\delta \hat{a}_\kappa
\\

&=
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa

\end{align*}
$$

と書くことができました。それでは肝心のシュレーディンガー方程式を解くにはどうするのでしょうか。

シュレーディンガー方程式を書いてみると、多体の波動関数を$\Phi(\tau_1,\tau_2,\cdots,\tau_N)$として、

$$
\begin{align*}


\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa
\Phi(\tau_1,\tau_2,\cdots,\tau_N)

=
E
\Phi(\tau_1,\tau_2,\cdots,\tau_N)

\end{align*}
$$

を満たす、固有値$E$と、多体波動関数の固有関数＝スレーター行列式$\Phi(\tau_1,\tau_2,\cdots,\tau_N)$を求める、という方程式になっています。

第二量子化表示で考える場合は、先ほど見つけた関係式：

$$
\begin{align*}
\hat{a}_\delta^\dagger\hat{a}_\delta
|\cdots \varphi_\delta\cdots|

&=
|\cdots\varphi_\delta\cdots|
\end{align*}
$$

つまり、生成・消滅演算子のラベルと同じ1電子状態を含むスレーター行列式は演算子$\hat{a}_\delta^\dagger\hat{a}_\delta$の固有状態になっていることを利用して、ハミルトニアンの行列を

$$
\begin{align*}
H_{\eta\zeta}
&=
\int\phi_\eta(\tau)\hat{H}\phi_\zeta(\tau)d\tau

=
\delta_{\eta\zeta}\varepsilon_\zeta

\end{align*}
$$

と対角化できる関数系を見つけて、多体の波動関数をこの関数系で展開することにすると、第二量子化表示のハミルトニアンを

$$
\mathcal{H} = \sum_\zeta\varepsilon_\zeta \hat{b}_\zeta^\dagger\hat{b}_\zeta
$$

と書くことができ、この時の1電子状態を並べたスレーター行列式$|\phi_\eta\phi_\omega\cdots\phi_\rho|$は、

$$
 \sum_\zeta\varepsilon_\zeta \hat{b}_\zeta^\dagger\hat{b}_\zeta
|\phi_\eta\phi_\omega\cdots\phi_\rho|
=
\left(
    \varepsilon_\eta + \varepsilon_\omega + \cdots \varepsilon_\rho
    \right)
|\phi_\eta\phi_\omega\cdots\phi_\rho|
$$

と、ハミルトニアン$\mathcal{H} = \sum_\zeta\varepsilon_\zeta \hat{b}_\zeta^\dagger\hat{b}_\zeta$の、固有値$E = (\varepsilon_\eta + \varepsilon_\omega + \cdots \varepsilon_\rho)$の固有状態になっていることがわかります。

というわけで、第二量子化表示で書き直したシュレーディンガー方程式を解くには、ハミルトニアンを

$$
\begin{align*}
H_{\eta\zeta}
&=
\sum_{\eta\zeta}\int\phi_\eta(\tau)\hat{H}\phi_\zeta(\tau)d\tau

=
\delta_{\eta\zeta}\varepsilon_\zeta

\end{align*}
$$

と対角化する関数系を求めて、それをスレーター行列式に詰めていきましょう。

…って、これじゃ結局今までと同じで、一体のハミルトニアンの固有状態を求めているだけじゃねえかと思うわけですが、まさにそうで、「第二量子化」等と仰々しい表現を使っていますが、本質的には今まで考えてきた「量子力学」と変わるところはなく、ただ表現が変わっただけなのでした。

とはいえ、ちょっとしたアプローチの違いはあります。
これまでは「波動関数を平面波で展開して係数を求める」とか、「波動関数をパラメータを含んだ試行関数で仮定し、変分法で近似的な波動関数を求める」というようにして波動関数の側に重きを置いたアプローチだったのですが、第二量子化表示の場合は先ほどの手順「ハミルトニアンを変形していき、自明な固有状態を持つ形（対角化された形）に変形する」のように、波動関数が議論の俎上に出てくることはあまりなくなり、主役が第二量子化表示されたハミルトニアンの方に移る印象です。


## 生成・消滅演算子の変換

最後に、ある関数系を生成・消滅する演算子から、別の関数系の演算子への変換規則について扱います。

細かい話は省略しますが、ある正規直交基底をなす関数系$\{\varphi_\nu\}$でハミルトニアンを

$$
\begin{align*}
\mathcal{H}
&=
\sum_{\delta\kappa}\int\varphi_\delta(\tau)\hat{H}\varphi_\kappa(\tau)d\tau\hat{a}^\dagger_\delta \hat{a}_\kappa
\\

&=
\sum_{\delta\kappa}H_{\delta\kappa}\hat{a}^\dagger_\delta \hat{a}_\kappa

\end{align*}
$$

と書いた後に、別の正規直交基底$\{\phi_\eta\}$での表現に書き換えたいことが良くあります。（具体例はTight-bindingモデルの第二量子化表示を考える際に、すぐに見ることになります）

これは、前節まで見てきたように多体の波動関数を完全正規直交$\{\varphi_\nu\}$で展開したという前提から、やっぱり別の完全正規直交系$\{\phi_\eta\}$で展開したと考えた場合、それらの関数の生成・消滅演算子の間にはどのような関係（変換規則）があるかを考えていくものです。

愚直には、（「シュレーディンガー方程式を解く」の節で考えたように）ゼロからハミルトニアンの行列を考えなおして新しい基底で$\int\phi_\eta(\tau)\hat{H}\phi_\zeta(\tau)d\tau$を計算しなおせばよいのですが、そうせずとも変換前の関数と変換後の関数の関係を使って、生成・消滅演算子を書き換えることができます。これは具体的な問題を解く際に良く使うことになる関係です。

### 基底の間の関係

ここで変更前の関数系と、変更後の関数系に対して、係数$c_{\mu\eta}$を用いて

$$
\varphi_\nu = \sum_\eta c_{\eta\nu}\phi_\eta
$$

また逆変換が

$$
\phi_\eta = \sum_\nu \tilde{c}_{\nu\eta}\phi_\eta
$$


と関係づけられているとします。なお、仮定よりどちらも正規直交基底をなすので、係数$c_{\mu\eta}$は

$$
c_{\eta\nu} = \int \phi_\eta^*(\tau)\varphi_\nu(\tau)d\tau\equiv
\braket{ \phi_\eta|\varphi_\nu}
\\
\tilde{c}_{\nu\eta} = \int \varphi_\nu^*(\tau)\phi_\eta(\tau)d\tau
=

c_{\eta\nu}^* = 
\braket{ \varphi_\nu|\phi_\eta}


$$

を満たし、この係数を上記のように$\braket{ \phi_\eta|\varphi_\nu}$等と定義します。（これはブラケット記法ですが、詳しいことはさておきここではとりあえず積分の略記法というくらいで使います）

また

$$
\varphi_\nu = \sum_\eta c_{\eta\nu}\phi_\eta
=
 \sum_\eta c_{\eta\nu}\sum_\zeta\tilde{c}_{\zeta\eta}\varphi_\zeta
$$

より、

$$
\sum_\eta\sum_\zeta c_{\eta\nu}\tilde{c}_{\zeta\eta} = \sum_\eta \sum_\zeta c_{\eta\nu} c_{\eta\zeta}^* = \sum_\zeta\delta_{\nu\zeta}
$$

または

$$
\sum_\eta\sum_\zeta\braket{ \phi_\eta|\varphi_\nu}\braket{ \varphi_\zeta|\phi_\eta}
=
 \sum_\zeta\delta_{\nu\zeta}

$$

または上式を並べ替えて、

$$
\sum_\eta \ket{\phi_\eta}\bra{\phi_\eta} = 1
$$

を満たします。（シンプルな例としてFourier級数展開$\varphi_\nu(\boldsymbol{r}) = \sum_{\boldsymbol{k}} c_{\boldsymbol{k}}^{\nu}e^{i\boldsymbol{k}\cdot\boldsymbol{r}}$や、Block関数からWannier関数への変換などが挙げられますが、これに限らず色々な変換を考えることになります。）

### 生成・消滅演算子の変換

さて、このような関係がある場合変更前の$\varphi_\nu$についての生成・消滅演算子と、変更後の$\phi_\eta$の生成・消滅演算子の間にはどのような対応があるでしょうか。

$\varphi_\nu$に対応した演算子を$\hat{a}_\nu,\hat{a}^\dagger_\nu$、変更後の$\phi_\eta$に対応した演算子を$\hat{b}_\eta, \hat{b}_\eta^\dagger$と置いて、両者の関係を考えてみます。


これは、何もない状態を$|\>\>|$と書いて以下のような関係を見ればよく、（普通はケットベクトルを使って$\ket{0}$等と書かれますが、どこかの脚注でも書いたようにブラケット記法について書くのがしんどいのでちょっとサボって変な書き方をしています）への作用を考えてみればよく、生成演算子については

$$
\hat{a}^\dagger_\nu|\>\>| = |\varphi_\nu| = \sum_\eta \braket{ \phi_\eta|\varphi_\nu}|\phi_\eta|=\sum_\eta \braket{ \phi_\eta|\varphi_\nu}\hat{b}_\eta^\dagger|\>\>|.
$$

より

$$
\hat{a}^\dagger_\nu
=
\sum_\eta c_{\eta\nu}b_\eta^\dagger
=
\sum_\eta \braket{ \phi_\eta|\varphi_\nu}\hat{b}_\eta^\dagger.
$$

また消滅演算子については

$$
|\varphi_\nu| = \sum_\eta \braket{ \phi_\eta|\varphi_\nu}|\phi_\eta|
$$

の対応から、それぞれを「消す」演算子を考えると、
左辺：

$$
\hat{a}_\nu|\varphi_\nu| =|\>\>|
$$

右辺：

$$
\begin{align*}
\sum_{\eta'}\braket{\varphi_\nu|\phi_\eta'}\hat{b}_{\eta'}\sum_\eta \braket{ \phi_\eta|\varphi_\nu}|\phi_\eta|

&=
\sum_{\eta'}\sum_\eta\braket{\varphi_\nu|\phi_\eta'} \braket{ \phi_\eta|\varphi_\nu}\hat{b}_{\eta'}|\phi_\eta|\\

&=
\sum_{\eta'}\braket{\varphi_\nu|\phi_\eta'}\sum_\eta \braket{ \phi_\eta|\varphi_\nu}\delta_{\eta\eta'}|\>\>|\\

&=
\sum_\eta\braket{\varphi_\nu|\phi_\eta} \braket{ \phi_\eta|\varphi_\nu}|\>\>|\\

&=|\>\>|

\end{align*}
$$

となるので、消滅演算子の変換は

$$
\hat{a}_\nu = \sum_\eta c_{\eta,\nu}^* \hat{b}_\eta = \sum_\eta \braket{\phi_\eta|\varphi_\nu} \hat{b}_\eta
$$

とすれば良さそうです。


以上まとめると、関数の変換

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

として変換できることがわかりました。この関係は次章のTight-bindingモデルの第二量子化表示について考える際に利用します。


## 1体演算子のその先へ

一旦、生成・消滅演算子演算子についての解説は以上です。

先ほど見たように1体演算子について扱う限りでは、結局やることは今までと全く変わらずただ書き直しただけ、という感じでした。
第二量子化表示が真価を発揮するのは2体演算子を扱ってからで、簡潔な表現と、様々な近似手法によって多体電子が、お互いに相互作用しあう複雑な状況を解析するための強力な道具になってくれます。

ただ、本稿ではあくまでTight-bindingモデルの第二量子化表示が最終目標ですので、2体演算子についてはさしあたり扱わないことにして、いずれ心に余裕ができたらまた色々と書いていこうと思います。

なお、2体演算子まで考えると、もはや波動関数について言及されることは無くなり、ひたすら生成・消滅演算子をファインマン図形で展開したり、生成・消滅演算子の平均値を求めたりすることで様々な物理量を求めていくことになる印象です。そのあたりは（量子）統計力学の分野になってきます。そこもいずれは書いていければなーなどと思っています。


## （余談）なんでこんなことするのか？

以上のように、多体の電子状態を扱うための第二量子化表示について簡単ですが学んできました。「本番」である二体演算子について言及せずできていないのは片手落ち感があるのですが、最後に余談としてなんでこんなことを考えるのか？についてだらだらと書き連ねておこうと思います。

個人的な理解としては、従来の「第一量子化」で問題を解くのと、「第二量子化」で解くのとでは、何かコンピュータで仕事をしたいと思ったときにどのプログラミング言語を使うかの違いに似ていると思っています。（別にプログラミングに詳しいわけではないので適当なことを書いてるかもですが）

例えば言語としてCを使っても、Pythonを使っても、最終的には機械語に翻訳されて、CPUやメモリの操作を行うわけで、「技術的には」作成できるプログラムの可能性は同じなのでしょう。しかしどちらも表現方法や得意な計算に特徴があり、特に機械学習関連のプログラムを書きたい場合はPythonには豊富なライブラリがあることから、Pythonを選択することが多くなりそうです。

また上記と関連するかもしれませんが、同じプログラムを書くにしても言語によって簡潔に表現できたり、込み入った表現が必要になったりするかと思います。


これと同じで「第一量子化」も「第二量子化」も、同じ方程式を書き換えているだけなので本質的に計算できることは同じです。
しかし、第二量子化の手法を選ぶと、これまで数多くの先人たちが積み上げてきた多体電子について扱う様々な手法を利用することができます。

また生成・消滅演算子の演算規則さえ受け入れてしまえば、かなりややこしい計算でもある程度機械的に進めていくことができます。

さらに、詳しくは書きませんが、「Habbardモデル」で導入される「On-Site Coulomb相互作用」すなわち同じ格子点上の電子間のみに働くCoulomb相互作用について、第二量子化表示では

$$
U\sum_in_{i,\uparrow}n_{i,\downarrow}
$$

と簡潔に書かれるのに対し、第一量子化表示で同じような相互作用の表現をするのはとても大変そうで、ちょっとすぐにはわかりません。

このように第二量子化表示を選ぶことには、様々な既存手法を利用できる点と、簡潔に考えている物理的な状態を表現できることに利点があるのだろうと考えております。

# おわりに

以上、私なりに丁寧に多体電子系のハミルトニアンを生成・消滅演算子で表す第二量子化表示について書きました。

改めて本章で得られた結果をまとめると、まず多体電子系のシュレーディンガー方程式

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

として変換できることが分かりました。


以上で第二量子化表示についての章を終わります。サクッとまとめるつもりが、やっぱり長くなってしまいました。次回でいよいよ最終章です。