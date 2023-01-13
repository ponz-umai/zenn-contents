---
title: "水素様原子のシュレディンガー方程式の具体的な計算"
emoji: "🎯"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# はじめに

水素原子中の電子の性質について、概要を[こちらの章](https://zenn.dev/ponzumai/articles/tight-binding-model-hydrogen-atom)でまとめたが、本章ではそこでは省略した具体的な計算やその他補足事項について述べる。

一定の量になったので現状で公開するが、内容は随時追記していく予定。

# シュレディンガー方程式の変形
## 極座標表示

前章で述べたように、一つの原子核に束縛された電子のシュレディンガー方程式は以下のように書ける。

$$
\left(-\frac{\hbar^2}{2m} \nabla^2  -\frac{Ze}{4\pi\epsilon _0r}  \right) \varphi(\boldsymbol{r}) =\epsilon \varphi (\boldsymbol{r} )
$$

この方程式はポテンシャルが$r$のみに依存することに着目して、$r$を独立した変数として扱うために極座標表示へ座標変換することで解析的に解くことができる。シュレディンガー方程式の座標変換について以下にまとめる。

極座標表示でのシュレディンガー方程式へ書き直すために、極座標で微分演算子$\nabla^2$を書き直すと、

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
## 変数分離
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

に分解する。
$\theta , \phi$ 部分はさらに変数分離ができて、まず$Y(\theta ,\phi) = \Theta(\theta )\Phi (\phi )$とおいて

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

これらの微分方程式は解析的な解をもつことが知られていて、それぞれ$R(r)$は「ラゲールの陪多項式」、$\Theta (\theta )$は「ルジャンドルの陪多項式」というやつで表される。また$Y(\theta,\phi) = \Theta\Phi$は、（に規格化の係数をかけたもの）「球面調和関数」と呼ばれる関数となる。

「陪多項式」とかなにやら物騒な名前の関数が出てきて、この辺で考えるのをやめた人も多いかもしれない（過去の私のように）。ただ実はそこまでややこしいことはしておらず、すごく大雑把に言うと
- 微分方程式の解が級数$\Theta (\theta ) = f(\theta )\sum_n a_n \cos^n\theta$の形で表せると仮定する
- 実際に微分方程式に代入して、方程式やその他条件（規格化や直交性、境界条件等）を満たすような$f(\theta ), a_i$の形を求める

ということをしているだけですので、あまりビビらなくてもOKです。上記のように級数解を仮定していることから、固有関数は多項式で表されることになる。

解の求め方は例えば[この本](https://www.amazon.co.jp/dp/4781910068)や、[このWEBページ](https://batapara.com/archives/legendre-differential-equation.html/)が参考になる。

以下でも解法を記載していく。（予定）

# 微分方程式の解法

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

## 角度部分
（追記予定）

## 動径部分
（追記予定）