---
title: "水素様原子のシュレディンガー方程式の具体的な計算"
emoji: "🙌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: []
published: false
---
todo:
水素原子のシュレディンガー方程式の解法については本記事では概要を日本語で書いて、具体的な式の形を描くことにして、具体的な解放はどっかで別記事でまとめよう。まずは今書いてる極座標の計算とか細かい話をこっちに移行する

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

となる。これは変数変換$\left( x,y,z \right)\rightarrow \left( r,\theta ,\phi \right)$ をした際の微分演算子の変換規則

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

を利用して30分くらい頑張ってA4用紙2枚分位の計算をすれば導出できる。なお、教科書によっては$\frac{\partial ^2}{\partial r^2}+ \frac{2}{ r} \frac{\partial }{\partial r}$の部分を$\frac{1}{r^2}\frac{\partial }{\partial r}\left( r^2\frac{\partial }{\partial r}  \right)$と書いてあることがあるが、これは

$$
\frac{1}{r^2}\frac{\partial }{\partial r}\left( r^2\frac{\partial }{\partial r}  \right) =
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
$$

であり同値となる。
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

となる。ここに$-\frac{2m}{\hbar^2} \frac{r^2}{RY}$をかけて整理して、

$$
\frac{r^2}{R}\left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        \right)
        +
        \frac{2m}{\hbar^2}r^2\left( \epsilon + \frac{e}{4\pi\epsilon _0r} \right)  
 =-\frac{\hat{\Lambda}Y }{Y} 
$$

と変数分離の形を作ることができる。両辺を定数$\lambda$とおいて、$r$に関する微分方程式

$$
-\frac{\hbar^2}{2m} \left(
        \frac{d ^2R}{d r^2} + \frac{2}{r}\frac{d R}{d r}
        -\frac{\lambda }{r^2} R
        \right)
        -
       \frac{e}{4\pi\epsilon _0r}R = \epsilon R
$$

と$\theta , \phi$ に関する微分方程式

$$
\hat{\Lambda }Y + \lambda Y = 0
$$

に分解することができます。
$\theta , \phi$ 部分はさらに変数分離ができます。まず$Y(\theta ,\phi) = \Theta(\theta )\Phi (\phi )$とおいて

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

両辺に$\sin^2\theta \frac{1}{\Theta \Phi }$をかけて整理して、

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

を得る。以上から、$r,\theta ,\phi$に関する微分方程式を満たす関数$R(r), \Theta (\theta ), \Phi (\phi )$と、その際の$\epsilon$を求めることができれば、水素原子中の電子のハミルトニアンに対するエネルギー定常状態の解$\varphi (\boldsymbol{r})=R(r)\Theta (\theta )\Phi (\phi )$と、その際の固有エネルギーを知ることができる。

これらの微分方程式は解析的な解をもつことが知られていて、それぞれ$R(r)$は「ラゲールの陪多項式」、$\Theta (\theta )$は「ルジャンドルの陪多項式」というやつで表される。また$Y = \Theta\Phi$は、「球面調和関数」と呼ばれる関数となる。

「陪多項式」とかなにやら物騒な名前の関数が出てきて、この辺で考えるのをやめた人も多いかもしれない（過去の私のように）。ただ実はそこまでややこしいことはしておらず、すごく大雑把に言うと
- 微分方程式の解が級数$\Theta (\theta ) = f(\theta )\sum_n a_n \cos^n\theta$の形で表せると仮定する
- 実際に微分方程式に代入して、方程式やその他条件（規格化や直交性、境界条件等）を満たすような$f(\theta ), a_i$の形を求める

ということをしているだけですので、あまりビビらなくてもOKです。上記のように級数解を仮定していることから、固有関数は多項式で表されることになるのですが、項の数は有限個にとどまります。

解の求め方は例えば[この本](https://www.amazon.co.jp/dp/4781910068)や、[このWEBページ](https://batapara.com/archives/legendre-differential-equation.html/)等を参照していただくとして、概要だけを以下に示すことにする。

改めて解くべき微分方程式を並べて書くと

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
まず、変数分離の過程で置いた定数$\lambda ,m$について以下の条件が課される：

$$
\lambda =l(l+1),\>\>\>\> l=0,1,2\dots,\\
\left|m\right|\leq l
$$
またこれらを用いて、$\Theta ,\Phi$は

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

と書かれる。$Y_l^m(\theta ,\phi )$は球面調和関数と呼ばれる関数系で、$l, m$に関して正規直交基底を張る。$P_l^m$はルジャンドル陪多項式と呼ばれる。

**落ち着いてほしい。**

多分これだけを見ても意味が分からないと思う。「ただの多項式です」とは何だったのか。そう思う気持ちもわかる。ただ上の式は一般項をスッキリ書くためにクソややこしくなっているだけなので落ち着いてほしい。まず$Y_l^m$の$(-1)^{\frac{|m|+m}{2} }\left[\frac{2l+1}{4\pi}\frac{(l-|m|)!}{(l+|m|)!} \right]^{1/2}$の部分はただの定数で、規格化のための係数と思えば怖くない。

また$e^{im\phi }$は微分方程式$\frac{d^2\Phi }{d\phi ^2} + m^2\Phi =0$の解であり$Y_l^m=\Theta \Phi$の$\Phi$部分に他ならない。というわけで残るは$P_l^m(\cos\theta )$部分なわけですが、これも$l-|m|$次の多項式を何やら微分しまくって一般項で表しているだけで、主要な部分のみとりだすと、まず二項定理で展開して

$$
(\omega ^2-1)^l = \omega ^{2l} + {}_l{\rm C}_{l-1}\omega ^{2l-2}(-1) + {}_l{\rm C}_{l-2}\omega ^{2l-4}(-1)^2 + \cdots
$$

を、

$$
\frac{d^l}{d\omega ^l}
$$

で$l$回微分して、第一項は$(l$の項$)\times \omega ^{2l-l}$、第二項は$(l$の項$)\times\omega ^{2l-1-l}\cdots$で、$\omega ^{2l-l-l}$の項以降は微分でゼロになる。というわけで最初の$l$回微分で「$l$の項」とした部分を$a_i(l)$とおくと

$$
a_1\omega ^l + a_2\omega ^{l-1} + \cdots a_l\omega ^0 
$$

というただの$\omega$の$l+1$項からなる多項式になるわけで、
その後

$$
\frac{d^{|m|}}{d\omega^{|m|}}
$$

で$|m|$回微分すると$\omega$の係数が$|m|$未満の項が無くなり最終的に$l+1-|m|$項の$\omega^{l-|m|},\omega ^{l-1-|m|}\cdots$の多項式に
$(1-\omega ^2)^{|m|/2}$掛けた関数になる、という構造になっておりまして、そう考えるとそこまでややこしい関数じゃないように思えてくるのではないでしょうか。

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

エネルギーはラベル$n$のみに依存し、$l$の値によらないことが特徴的だが、これはクーロンポテンシャルの場合のみで、一般の球対称ポテンシャルの場合は$l$による。（例えばこの後扱う多電子原子だとポテンシャルの形がクーロンポテンシャルからずれ、その結果エネルギーが$l$にも依存するようになる。）

関数$R_{nl}$は、通常$l$の部分を先ほど述べた$s,p,d\cdots$と置き換えて$R_{1s}, R_{2p}\cdots$などと表現される。

以上まとめて、エネルギー固有状態を表す波動関数はラベル$n,l,m$によって指定される関数$\varphi_{nlm}=R_{nl}(r)Y_l^m(\theta ,\phi )$となる。

例えば、波動関数$R_n$において $n=1$の場合は$l$の取りうる値は
一般項をはじめに書いても意味が分からないと思うので最初の数項のみ具体的に以下に示すと、

$$
R_{}
$$

これらの関数はそれぞれ正規化され、異なる関数同士は直交している。また全てのラベルについて集めると完全系をなしており、正規直交基底となっている。