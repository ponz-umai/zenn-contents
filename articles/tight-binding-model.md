---
title: "tight binding入門～モデルの成り立ちを理解する"
emoji: "📚"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: []
published: false
---
# はじめに
## ざっくりtight bindingの説明
物性物理では以下のようなtight-binding modelがよく使われます。
$$
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }\left( -t_{ij}a^\dagger_{i\sigma }a_{j\sigma } + \mathrm{h.c}. \right) 
$$
ここで$ \left<i,j \right>$は最近接の格子点について和を取ることを表します。$t_{ij}$は飛び移り積分、$a^\dagger_{i\sigma },a_{j\sigma }$はそれぞれ格子点$i(j)$、スピン$\sigma $の電子の生成（消滅）演算子で、$\rm{h.c.}$は第一項のエルミート共役です。

このモデルは固体内の電子が格子点（$\simeq $原子核）に強く束縛され、局在していると仮定した描像に対応しています。

tight-binding modelは非常にシンプルな構造ですが、そこに電子間のオンサイトクーロン相互作用を加えたHubbard model
$$
\mathcal{H} = \sum_{\left <i,j \right>,\sigma }\left( -t_{ij}a^\dagger_{i\sigma }a_{j\sigma } + \mathrm{h.c}. \right)  + U\sum_in_{i\uparrow}n_{i\downarrow}
$$
に始まり、局在磁気モーメントを加えたAnderson model等々様々なモデルの出発点となっており、それらからMot絶縁体やBCS超伝導、Anderson局在、はたまた今ホットなトポロジカル物性、量子ビット候補のマヨラナ粒子、、等々様々な物理現象が描き出されることになり、物性物理の勉強・研究をしていれば必ず出会うことになります。または物性とは異なる分野の方でも目にする機会はあるのではないかと思います。

## っていうことなのですが
ひとたび上記のような生成・消滅演算子で記述されたモデル（ハミルトニアン）を受け入れてしまえば、その先は演算子の交換関係や、平均場近似等々を使ってある程度機械的にモデルの解析を進めていくことができます。

しかし「サイト$i,j$の電子って何？」とか、「飛び移り積分って？」とか「電子を作ったり消したりする演算子っていったい何？？？」とか考えてみるとなんだかよくわからない、（学部時代を適当に過ごしてきた）大学院生や、他分野の研究者・技術者の方もいるのではないでしょうか。
というか何を隠そう、私が学部時代に量子力学を適当に勉強した状態で修士に進み、その辺にもやもやを抱えながら研究してきた落ちこぼれ院生だったのですが（そしてそれを解消できないまま修了しちゃったのですが）、このままだと死んでも死にきれないので過去の自分がブラックボックスにしてきた部分を理解するべく勉強した内容をまとめたノートがこちらです。

量子力学についてのなんとなくの理解を前提に、上記tight-binding modelの第二量子化表示について納得できるまでの内容をまとめます。とはいえ全部書くとすごい量になってしまい終わらなさそうなので、適宜既存の教科書に説明を譲りつつ、教科書では省略されていたりあまり触れられていない部分を中心にまとめていきます。質問やわかりにくい点などあればコメントを頂けますと幸いです。（もちろん単なる感想でもいただけるととても喜びます。）

また、上述の通り落ちこぼれ修士が独学をベースに書いているノートなので、勘違いや間違いは多数含まれている可能性がありますし、おかしな書き方をしている部分もあるかと思いますのでご注意ください（指摘していただければ嬉しいです）。

# 量子力学の復習
tight-bindingモデルは、下図のように、「固体は多数の原子核と電子が周期的に配列している」「電子はほぼ原子核の近くに局在しており（原子核に束縛されており）、固体を構成する原子（または分子）が孤立している場合と近い状態にある」という状況を想定したモデルです。

したがってその出発点としてまず、孤立した一原子（原子核＋（複数の）電子中の、電子の振る舞いについて理解することとなります。
一つの原子核＋一つの電子すなわち「水素様原子」の電子の状態は厳密に求めることができますが、電子が二つ以上（多電子系）になると近似が必要となり、また、多電子系特有の考え方が必要になってきます。

以下では、そもそも電子の振る舞いを記述するのに必要な量子力学の考え方（シュレディンガー方程式）から出発し、固体構成要素となる1原子中の1電子の振る舞い、そして多電子になった場合の扱いについて簡単にまとめます。

## シュレディンガー方程式
本稿で主に扱う電子や、ちらっと出てくる原子核（陽子と中性子）、光（光子）等ミクロな粒子の振る舞いは、以下のようなシュレディンガー方程式によって記述される。
$$
i\hbar \frac{\partial}{\partial t} \psi\left( {\boldsymbol r}, t\right) 
= \hat{H}\left( \hat{\boldsymbol{r}},\hat{\boldsymbol{p}} \right) \psi\left( \boldsymbol{r} ,t\right)
$$
ここで$ \hat{H}\left( \hat{\boldsymbol{r}},\hat{\boldsymbol{p}} \right) $は系のハミルトニアンに対して、$ \hat{\boldsymbol{r}} \rightarrow \boldsymbol{r},\hat{\boldsymbol{p}} \rightarrow -i\hbar \nabla  $と置き換えたもので、$　\psi\left( \boldsymbol{r} ,t\right)$に対して作用します。
$　\psi\left( \boldsymbol{r} ,t\right)$は波動関数と呼ばれ、絶対値の二乗$\left|\psi\left( \boldsymbol{r} ,t\right)\right|^2$が、「時刻$t$, 位置$\boldsymbol r$で粒子が観測される確率」（正確には確率密度）に対応しています。今のところは$　\psi\left( \boldsymbol{r} ,t\right)$は一粒子だけ考えることにして、対応してハミルトニアン演算子$\hat{H}$も一粒子のハミルトニアンとします。（本稿では一粒子に作用する演算子を$\hat{H}$のように斜体にハットをつけた文字で表し、多粒子の演算子を$\mathcal{H}$のように筆記体で表すことにします）

時間によって変化しない状態を求める場合は$\psi\left( \boldsymbol{r} ,t\right) = \varphi\left( \boldsymbol{r} \right)e^{-i\omega t} $と置いてシュレディンガー方程式に代入し、
$$
i\hbar \varphi(r) (-i\omega) e^{-i\omega t} = e^{-i\omega t}\hat{H}\varphi(r) \\
\Rightarrow \hbar \omega = \hat{H}\varphi(r)　(両辺をe^{-i\omega t}で割った)
$$
ここで$\hbar \omega$ はドブロイの関係式より、物質のエネルギーと解釈できるので、これを$\epsilon $と置いて「時間に依存しないシュレディンガー方程式」

$$
\hat{H}\varphi\left( \boldsymbol{r} \right) = \epsilon \varphi\left( \boldsymbol{r} \right)
$$
にたどり着きます。これは演算子$\hat{H}$の固有値方程式となっており、$\hat{H}$の固有関数を求めることにより、ハミルトニアンによって設定された状況に対して実現する、粒子の（時間変化しない確率密度の）分布を関数として得ることができます。

以下では基本的に時間変化しない解を考えることにして、上記「時間に依存しないシュレディンガー方程式」を単に「シュレディンガー方程式と書きます。

## 水素原子中の電子のハミルトニアンと波動関数
先述の通り固体は、無数の原子が集まった状態と考えられ、ハミルトニアン（系の全エネルギー）は

（各電子の運動エネルギー）＋　（各原子核の運動エネルギー）＋（電子と原子核の間の（引力）相互作用ポテンシャル）＋（原子核間の（斥力）相互作用ポテンシャル）＋（電子間の相互作用ポテンシャル）＋（その他色々）

のような要素で構成されるはずです。それら全てを取り扱うのは辛すぎるので様々な近似をしていくことになるのですが、まずは最も簡単なパターンとして原子核1個＋電子1個の水素原子のシュレディンガー方程式を考えます。これは近似無しで解くことができて、固体＝多電子・多原子の状態を考えていくための出発点となります。

原子核の電荷を$Ze$（水素原子の場合は$Z=1$）と電子間のクーロン相互作用が電子の位置を$\boldsymbol{r}$, 原子核の座標を$\boldsymbol{R}$として、
$$ V({\boldsymbol r})=-\frac{Ze^2}{4\pi \epsilon_0 \left|{\boldsymbol r - R}\right|}
$$
原子核の位置$\boldsymbol{R}$を原点に選ぶと、ポテンシャルが$\boldsymbol{r}$の絶対値$r\geq  0)$のみの関数となることを踏まえて、シュレディンガー方程式は
$$
\left(-\frac{\hbar^2}{2m} \nabla^2  -\frac{Ze}{4\pi\epsilon _0r}  \right) \varphi(\boldsymbol{r}) =\epsilon \varphi (\boldsymbol{r} )
$$
となります。
### 水素原子中の電子の固有状態
この微分方程式は解析的に解けて、先にその解について記載します。

まず固有関数は3つのラベル（量子数）$n, l, m$により指定される動径$r$部分の関数$R_{nl}(r)$、角度$\theta $部分の関数$P_l^m(\theta )$（に規格直交化の係数がかかったもの）、角度$\phi $部分の関数$\frac{1}{\sqrt{2\pi}}e^{im\phi}$ の積で表され、ラベル$n,l,m$の間には以下の関係があります。
$$
n=1,2,3\cdots,\\
l: 0または正の整数, かつ \>\>\>\> l+1\leq n\\
m: 整数、かつ\>\>\>\> |m|<l
$$
角度$\theta$部分の関数と$\phi $部分の関数を、規格化等の係数も全て合わせて$Y_l^m(\theta ,\phi )$と書き、この関数は規格直交します。

固有エネルギーは量子数$n$のみに依存し、

$$
\epsilon _n = -\frac{me^4}{(4\pi\epsilon )^22\hbar^2} \frac{1}{n^2} , \>\>\>\> n=1,2,\cdots
$$
となります。ある$n$を選ぶと、その値に対して$l$が$l=0,1,\cdots n-1$の$n$個あり、かつそれぞれの$l$に対して$m$が$m=-l, -(l-1)\cdots ,0,\cdots l$の$2l+1$個あるので、合計$\sum_{k=0}^{n-1}(2k+1)=n^2$個の縮退した状態が存在することになります。

この事情はポテンシャルエネルギーがクーロンエネルギーの形をしていることによる特殊なもので、ポテンシャルがクーロンエネルギーからずれる場合はエネルギーが$l$にも依存するようになります。（とはいえずれが小さい場合は大体同じ程度のエネルギーにとどまると考えられます。この性質はTight-binding近似を考える際に出てきますので、そこでまた思い出していただければと思います。）

一般項は無駄にややこしいので置いておいて先に具体的な解の形を書くと、
#### n=1
$$
R_{10}(r) = \left( \frac{Z}{a_0}  \right)^{3/2}2e^{-Zr/a_0}\\
Y_0^0 = 
$$

#### n=2
$$
R_{20}(r) = \left( \frac{Z}{a_0}  \right)^{3/2}\frac{1}{\sqrt{2}}\left( 1-\frac{1}{2} \frac{Zr}{a_0}  \right)  e^{-Zr/2a_0}\\


$$


### 極座標表示のシュレディンガー方程式
この方程式はポテンシャルが$r$のみに依存することに着目して、$r$を独立した変数として扱うために極座標表示へ座標変換することで解析的に解くことができる。極座標表示でのシュレディンガー方程式へ書き直すために、極座標で微分演算子$\nabla^2$を書き直すと、
$$
\nabla^2 = \frac{\partial ^2}{\partial r^2} + \frac{2}{r}\frac{\partial }{\partial r} + \frac{1}{r^2}\hat{\Lambda} (\theta ,\phi), \\
\hat{\Lambda} (\theta ,\phi) = \left\{
     \frac{1}{\sin \theta}  \frac{\partial }{\partial \theta } \left( \sin \theta \frac{\partial }{\partial \theta }  \right) 
     +
     \frac{1}{\sin ^2\theta }\frac{\partial ^2}{\partial \phi ^2}  
\right\}     
$$
となる。これは変数変換$\left( x,y,z \right)\rightarrow \left( r,\theta ,\phi \right) $ をした際の微分演算子の変換規則
$$
\frac{\partial }{\partial x}  = \frac{\partial r}{\partial x}\frac{\partial }{\partial r}  + 
\frac{\partial \theta }{\partial x}\frac{\partial  }{\partial \theta} 
+
\frac{\partial \phi}{\partial x}\frac{\partial }{\partial \phi}
$$
（$y,z$についても同様）
と、$x,y,z$と$r,\theta ,\phi$の間の関係式
$$
r = x^2+y^2+z^2,\\
\tan^2\theta = \frac{x^2+y^2}{z^2} ,\\
\tan\phi=\frac{y}{x} 
$$
を利用して30分くらい頑張ってA4用紙2枚分位の計算をすれば導出できる。なお、教科書によっては$\frac{\partial ^2}{\partial r^2}+ \frac{2}{ r} \frac{\partial }{\partial r} $の部分を$\frac{1}{r^2}\frac{\partial }{\partial r}\left( r^2\frac{\partial }{\partial r}  \right)$と書いてあることがあるが、これは
$\frac{1}{r^2}\frac{\partial }{\partial r}\left( r^2\frac{\partial }{\partial r}  \right) =
 \frac{1}{r^2}
\left[
    \left(
        \frac{\partial }{\partial r}r^2
    \right)\frac{\partial }{\partial r} 
    +
    r^2
    \frac{\partial }{\partial r}\frac{\partial }{\partial r} 
    \right]
    =\frac{\partial ^2}{\partial r^2}+ \frac{2}{ r} \frac{\partial }{\partial r} 
     $であり同値となる。
移動中などの隙間時間に一度やってみると納得できると思う。

上式をハミルトニアンに代入して、シュレディンガー方程式は
$$
\left\{
    -\frac{\hbar^2}{2m}\left(
        \frac{\partial ^2}{\partial r^2} + \frac{2}{r}\frac{\partial }{\partial r} + \frac{1}{r^2}\hat{\Lambda} (\theta ,\phi)
         \right) -\frac{e}{4\pi\epsilon _0r}  
    \right\}
 \varphi(\boldsymbol{r}) =\epsilon \varphi (\boldsymbol{r} ),\\
\hat{\Lambda} (\theta ,\phi) = \left\{
     \frac{1}{\sin \theta}  \frac{\partial }{\partial \theta } \left( \sin \theta \frac{\partial }{\partial \theta }  \right) 
     +
     \frac{1}{\sin ^2\theta }\frac{\partial ^2}{\partial \phi ^2}  
\right\} 
$$
となる。
### 変数分離
両辺に$r^2$をかけることで「$r$のみに依存する部分」と「$\theta ,\phi$のみに依存する部分」に分けられそうな形になった。そこで$\varphi(\boldsymbol{r} ) = R(r)Y(\theta ,\phi)$とおいて代入し、
$$
-\frac{\hbar^2}{2m}\left\{
    \left(
        \frac{\partial ^2R}{\partial r^2} + \frac{2}{r}\frac{\partial R}{\partial r}
        \right)Y
         + \frac{R}{r^2}\hat{\Lambda} Y
    \right\}      
         -\frac{e}{4\pi\epsilon _0r}  RY
 =\epsilon RY
$$
となる。ここに$-\frac{2m}{\hbar^2} \frac{r^2}{RY} $をかけて整理して、
$$
\frac{r^2}{R}\left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        \right)
        +
        \frac{2m}{\hbar^2}r^2\left( \epsilon + \frac{e}{4\pi\epsilon _0r} \right)  
 =-\frac{\hat{\Lambda}Y }{Y} 
$$
と変数分離の形を作ることができる。両辺を定数$\lambda $とおいて、$r$に関する微分方程式
$$
-\frac{\hbar^2}{2m} \left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        -\frac{\lambda }{r^2} R
        \right)
        -
       \frac{e}{4\pi\epsilon _0r}R = \epsilon R
$$
と$\theta , \phi $ に関する微分方程式
$$
\hat{\Lambda }Y + \lambda Y = 0
$$
に分解することができます。
$\theta , \phi $ 部分はさらに変数分離ができます。まず$Y(\theta ,\phi) = \Theta(\theta )\Phi (\phi ) $とおいて
$$
\left\{
     \frac{1}{\sin \theta}  \frac{\partial }{\partial \theta } \left( \sin \theta \frac{\partial }{\partial \theta }  \right)\Theta \right\}
     \Phi  
     +
     \frac{1}{\sin ^2\theta }\Theta\frac{\partial ^2}{\partial \phi ^2}  
      \Phi 
      +
      \lambda \Theta \Phi=0 
$$
両辺に$\sin^2\theta \frac{1}{\Theta \Phi } $をかけて整理して、
$$
\frac{1}{\Theta } \left\{
     \sin \theta \frac{\partial }{\partial \theta } \left( \sin \theta \frac{\partial }{\partial \theta }  \right)\Theta
     \right\}
     +
     \lambda \sin^2\theta 
     =

    - \frac{1}{\Phi}\frac{\partial ^2}{\partial \phi ^2}  
      \Phi 
$$
両辺を定数$m^2$とおくと、
$$
\frac{1}{\sin\theta } \frac{d}{d\theta } \left( \sin\theta \frac{d\Theta }{d\theta }   \right) + \left( \lambda - \frac{m^2}{\sin^2\theta }  \right) \Theta =0,\\
\frac{d^2\Phi }{d\phi ^2} + m^2\Phi =0
$$
を得る。以上から、$r,\theta ,\phi $に関する微分方程式を満たす関数$R(r), \Theta (\theta ), \Phi (\phi )$と、その際の$\epsilon $を求めることができれば、水素原子中の電子のハミルトニアンに対するエネルギー定常状態の解$\varphi (\boldsymbol{r})=R(r)\Theta (\theta )\Phi (\phi )$と、その際の固有エネルギーを知ることができる。

これらの微分方程式は解析的な解をもつことが知られていて、それぞれ$R(r)$は「ラゲールの陪多項式」、$\Theta (\theta )$は「ルジャンドルの陪多項式」というやつで表される。また$Y = \Theta\Phi $は、「球面調和関数」と呼ばれる関数となる。

「陪多項式」とかなにやら物騒な名前の関数が出てきて、この辺で考えるのをやめた人も多いかもしれない（過去の私のように）。ただ実はそこまでややこしいことはしておらず、すごく大雑把に言うと
- 微分方程式の解が級数$\Theta (\theta ) = f(\theta )\sum_n a_n \cos^n\theta $の形で表せると仮定する
- 実際に微分方程式に代入して、方程式やその他条件（規格化や直交性、境界条件等）を満たすような$f(\theta ), a_i$の形を求める

ということをしているだけですので、あまりビビらなくてもOKです。上記のように級数解を仮定していることから、固有関数は多項式で表されることになるのですが、項の数は有限個にとどまります。

解の求め方は例えば[この本](https://www.amazon.co.jp/dp/4781910068)や、[このWEBページ](https://batapara.com/archives/legendre-differential-equation.html/)等を参照していただくとして、概要だけを以下に示すことにする。

改めて解くべき微分方程式を並べて書くと、
$$
-\frac{\hbar^2}{2m} \left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        -\frac{\lambda }{r^2} R
        \right)
        -
       \frac{e}{4\pi\epsilon _0r}R = \epsilon R,\\
\frac{1}{\sin\theta } \frac{d}{d\theta } \left( \sin\theta \frac{d\Theta }{d\theta }   \right) + \left( \lambda - \frac{m^2}{\sin^2\theta }  \right) \Theta =0,\\
\frac{d^2\Phi }{d\phi ^2} + m^2\Phi =0.
$$
### 角度部分のシュレディンガー方程式
まず、変数分離の過程で置いた定数$\lambda ,m $について以下の条件が課される：
$$
\lambda =l(l+1),\>\>\>\> l=0,1,2\dots,\\
\left|m\right|\leq l
$$
またこれらを用いて、$\Theta ,\Phi $は
$$
\Theta(\theta ) \Phi(\phi )= Y_l^m(\theta ,\phi ),\\
Y_l^m(\theta ,\phi )= 
    (-1)^{\frac{|m|+m}{2} }
        \left[
            \frac{2l+1}{4\pi}\frac{(l-|m|)!}{(l+|m|)!}  
        \right]^{1/2}P_l^m(\cos\theta )e^{im\phi },\\
P_l^m(\omega ) = (1-\omega ^2)^{|m|/2}\frac{d^{|m|}}{d\omega^{|m|}} P_l(\omega ),\\
 P_l(\omega ) = \frac{1}{2^ll!}\frac{d^l}{d\omega ^l} \left( \omega ^2 -1 \right) ^l
$$
と書かれる。$ Y_l^m(\theta ,\phi )$は球面調和関数と呼ばれる関数系で、$l, m$に関して正規直交基底を張る。$P_l^m$はルジャンドル陪多項式と呼ばれる。

**落ち着いてほしい。**

多分これだけを見ても意味が分からないと思う。「ただの多項式です」とは何だったのか。そう思う気持ちもわかる。ただ上の式は一般項をスッキリ書くためにクソややこしくなっているだけなので落ち着いてほしい。まず$Y_l^m$の$(-1)^{\frac{|m|+m}{2} }
        \left[
            \frac{2l+1}{4\pi}\frac{(l-|m|)!}{(l+|m|)!}  
        \right]^{1/2}$の部分はただの定数で、規格化のための係数と思えば怖くない。

また$e^{im\phi } $は微分方程式$\frac{d^2\Phi }{d\phi ^2} + m^2\Phi =0$の解であり$Y_l^m=\Theta \Phi $の$\Phi $部分に他ならない。というわけで残るは$P_l^m(\cos\theta )$部分なわけですが、これも$l-|m|$次の多項式を何やら微分しまくって一般項で表しているだけで、主要な部分のみとりだすと、まず二項定理で展開して
$$
(\omega ^2-1)^l = \omega ^{2l} + {}_l{\rm C}_{l-1}\omega ^{2l-2}(-1) + {}_l{\rm C}_{l-2}\omega ^{2l-4}(-1)^2 + \cdots
$$
を、
$$
\frac{d^l}{d\omega ^l}
$$
で$l$回微分して、第一項は$(l$の項$)\times \omega ^{2l-l}
$、第二項は$(l$の項$)\times\omega ^{2l-1-l}\cdots$で、$\omega ^{2l-l-l}$の項以降は微分でゼロになる。というわけで最初の$l$回微分で「$l$の項」とした部分を$a_i(l)$とおくと
$$
a_1\omega ^l + a_2\omega ^{l-1} + \cdots a_l\omega ^0 
$$
というただの$\omega $の$l+1$項からなる多項式になるわけで、
その後
$$
\frac{d^{|m|}}{d\omega^{|m|}}
$$
で$|m|$回微分すると$\omega $の係数が$|m|$未満の項が無くなり最終的に$l+1-|m|$項の$\omega^{l-|m|},\omega ^{l-1-|m|}\cdots $の多項式に$
(1-\omega ^2)^{|m|/2}
$掛けた関数になる、という構造になっておりまして、そう考えるとそこまでややこしい関数じゃないように思えてくるのではないでしょうか。

また、$l$と$m$の関係についても多項式の項数$l+1-|m|$が正の値を取らないと式も何もなくなってしまうので、$l\geq |m|$という条件が出てくるのも納得できるかと思います。
とはいえごちゃごちゃと抽象的なことを書いてきましたので、具体的な$l,m$について$Y_l^m$の関数系を以下に書いてみます。これを見ると本当にただの多項式だとわかっていただけるかと思います。

#### $l=0$の場合（$m$に許される値は$0$のみで、$Y_0^0$のみ）
$$
Y_0^0 = \frac{1}{\sqrt{4\pi}} 
$$

#### $l=1$の場合：$m=-1, 0,1$を取り得て、$Y_1^-1, Y_1^0,Y_1^1$を取る
$$
Y_1^0 = \sqrt{\frac{3}{4\pi} }\cos\theta ,\\

$$
- ・・・

というような関数系$Y_l^m$が得られる。$l$を「方位量子数」、$m$を「磁気量子数」と呼ぶ。また$l=0,1,2,3\cdots$の状態をそれぞれ$s,p,d,f\cdots$状態と呼ぶ。

### 動径部分

$R$に関する微分方程式の解はラベル$n$によって分類され、また（$\lambda =l(l+1)$と置いたことから）方位量子数$l$にも寄る固有関数$R_{nl}$と、$n$のみに依存する固有値（固有エネルギー）
$$
\epsilon _n = -\frac{me^4}{(4\pi\epsilon )^22\hbar^2} \frac{1}{n^2} , \>\>\>\> n=1,2,\cdots
$$
となる。
ここで微分方程式を解く過程で$n$と$l$の間に
$$
n\geq l+1
$$
の制約が出てきて、$n$に対する$l$の値の上限が設定されることになる。

エネルギーはラベル$n$のみに依存し、$l $の値によらないことが特徴的だが、これはクーロンポテンシャルの場合のみで、一般の球対称ポテンシャルの場合は$l $による。（例えばこの後扱う多電子原子だとポテンシャルの形がクーロンポテンシャルからずれ、その結果エネルギーが$l$にも依存するようになる。）

関数$R_{nl}$は、通常$l$の部分を先ほど述べた$s,p,d\cdots$と置き換えて$R_{1s}, R_{2p}\cdots$などと表現される。

以上まとめて、エネルギー固有状態を表す波動関数はラベル$n,l,m$によって指定される関数$\varphi_{nlm}=R_{nl}(r)Y_l^m(\theta ,\phi )$となる。

例えば、波動関数$R_n $において $n=1$の場合は$l$の取りうる値は
一般項をはじめに書いても意味が分からないと思うので最初の数項のみ具体的に以下に示すと、
$$
R_{}
$$
これらの関数はそれぞれ正規化され、異なる関数同士は直交している。また全てのラベルについて集めると完全系をなしており、正規直交基底となっている。



## 多電子原子中の電子
## 多電子状態の波動関数（スレーター行列式）

### Hartree Fock近似

# 固体（結晶）中の電子状態、そしてtight-binding model（第一量子化）
前章では原子一つに対する多電子の取り扱いについて整理しました。本章ではその内容を踏まえて、多数の原子が集まり、その結果多数の原子核と多数の電子が存在する固体内の電子状態について、tight-binding近似を行い、tight-bindingモデルの基本的な考え方・近似の入れ方を整理します。

## 結晶内の電子の性質
### ブロッホ関数・結晶運動量
### ブリルアンゾーンとエネルギーバンド

## 周期的な原子ポテンシャルの結晶ハミルトニアン
### 原子ポテンシャルの仮定
### Born-Oppenheimer近似
### 一電子ハミルトニアン
## tight-binding近似
### 近似①：ブロッホ関数の原子軌道関数による近似
### 変分法による基底状態の導出
### 近似②：異なるサイト間の重なりが小さいという近似
### 簡単な例とエネルギーバンド図

# tight-binding model の第二量子化
ここまでで「$t_{ij}$って何？」という疑問に対して大体のイメージがついたかと思います。また、「サイト$i,j$の電子って何？」という疑問についても、大体わかったと思います。

いよいよ最後に、「電子を作ったり消したりする演算子って何？」という疑問を解消するために、第二量子化表示のtight-bindingハミルトニアンについて本章でまとめていきます。
## 生成・消滅演算子を用いたシュレディンガー方程式の書き換え（第二量子化表示）

### 一体演算子の書き換え
### (参考)二体演算子の書き換え
### 固体内の電子のハミルトニアンの第二量子化表示
### 射影演算子による有効ハミルトニアンの導出

# おわりに








