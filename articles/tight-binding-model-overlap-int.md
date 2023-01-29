---
title: "重なり積分（Overlap Integral）の物理的意味（の妄想）"
emoji: "📚"
type: "idea" # tech: 技術記事 / idea: アイデア
topics: ["quantum","quantumcomputing","quantumcomputer","物理","物理学"]
published: true
---
# はじめに

本章では前章に続いて、Tight-bindingモデルに突然出てくる重なり積分：

$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r}
$$

の物理的意味について考えていきます。というか妄想していきます。ここで、この後に言及しやすくするためにこの積分を


$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r}

\equiv

s_{\boldsymbol{R},\boldsymbol{R}'}^{m,l}
$$

と定義しておきます。（重なり積分は大体無慈悲にゼロにされるので、あまりこのような記法は見かけませんが、重なり積分は大体記号$s$で表されるのでそのあたりを取り入れて本稿ではこのような記号を導入します）


特に、Tight-binding近似においては先述のように往々にして

$$
\int\varphi_m^*(\boldsymbol{r}-\boldsymbol{R}) \varphi_l(\boldsymbol{r-\boldsymbol{R}'})d\boldsymbol{r} \simeq \delta_{m,l}\delta_{\boldsymbol{R},\boldsymbol{R}'}
$$

と、異なる格子点間・軌道間でゼロと近似される不遇な量ですが、飛び移り積分と同じく「何を計算しているか」は分かるものの「どういう物理的なイメージに対応しているのか」が、なんとなくわかるようでわかりません。我々は一体何を「小さい」と考えているのでしょうか？

波動関数の自乗は確率分布を表す、というのは量子力学の最初の方に学ぶわけですが、異なる波動関数（一方は複素共役）の積の積分が一体何を表すのか、一応正しいと考えられていることをもとに頑張って考えてみることにしました。

# 重なり積分の物理的イメージ

さて、[前章](https://zenn.dev/ponzumai/articles/tight-binding-model-hopping-int)の飛び移り積分は、「孤立した原子軌道に結晶のハミルトニアンを作用させた際の展開係数」として得られ、それを波動関数の時間発展と結びつけることができたわけでした。

一方重なり積分はどのような物理的イメージと結びつく操作に対応して出てくるでしょうか？

## よく見る説明（量子化学などで）

まず、よく見る説明として「線形結合を取った波動関数

$$
\psi(\boldsymbol{r}) = C(\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}') + \phi_{m}(\boldsymbol{r}-\boldsymbol{R}))
$$

で表される電子の状態（$C$は規格化定数）があったとして、左側の波動関数$\phi_{m'}(\boldsymbol{r}-\boldsymbol{R}')$と右側の波動関数$\phi_{m}(\boldsymbol{r}-\boldsymbol{R})$が「重なっている」部分に電子がいる確率」のような物があります。（特に量子化学などで分子について言及される辺りで）

特に、分子軌道の導入として水素分子の波動関数に関して言及されることが多いので、そのような状況を念頭に置いて、孤立水素原子の$1s$原子軌道を$\chi(\boldsymbol{r})$と置き、水素原子A,Bを考えてそれぞれの水素原子の$1s$原子軌道を$\chi_A(\boldsymbol{r}), \chi_B(\boldsymbol{r})$と置きます。

水素原子$A,B$が分子を形成した際の電子の波動関数を

$$
\psi(\boldsymbol{r}) = C(\chi_A(\boldsymbol{r}) + \chi_B(\boldsymbol{r}))
$$

と近似的に表すと、この波動関数を測定した際の位置の分布を表す確率密度$\rho(\boldsymbol{r}) = |\psi(\boldsymbol{r}) |^2$は、

$$
\rho(\boldsymbol{r}) = |C|^2
\left(\chi_A^2(\boldsymbol{r}) + \chi_B^2(\boldsymbol{r}) + 2\chi_A(\boldsymbol{r}) \chi_B(\boldsymbol{r})
\right)
$$

です。ここで$1s$軌道関数が実関数であることを反映しています。

確率密度の全空間での積分$\int\rho(\boldsymbol{r})dr$は$1$なので、定積分の積分範囲を省略し、また係数$C$も実数とすると、

$$
1 = C^2\int
\left(
     \chi_A^2(\boldsymbol{r}) + |\chi_B^2(\boldsymbol{r}) + 2\chi_A(\boldsymbol{r}) \chi_B(\boldsymbol{r}) 
     \right)
     d\boldsymbol{r}
$$

ですが、ここで「重なり積分」

$$
S \equiv \int \chi_A(\boldsymbol{r}) \chi_B(\boldsymbol{r}) 
     d\boldsymbol{r}
$$

を定義すると、$\int \chi^2(\boldsymbol{r})d\boldsymbol{r}=1$より、$C$は正負どちらにとっても同じなので正にとり

$$
C = \frac{1}{\sqrt{2(1 + S)}}
$$

となり、確率密度は

$$
\begin{align*}
\rho(\boldsymbol{r}) &=
 \frac{1}{\sqrt{2(1 + S)}}
\left(\chi_A^2(\boldsymbol{r}) + \chi_B^2(\boldsymbol{r}) + 2\chi_A(\boldsymbol{r}) \chi_B(\boldsymbol{r})
\right)\\

&=
 \frac{1}{\sqrt{2(1 + S)}}
\left(\chi_A^2(\boldsymbol{r}) + \chi_B^2(\boldsymbol{r})
\right)

+
\sqrt{\frac{2}{1 + S}}
 \chi_A(\boldsymbol{r}) \chi_B(\boldsymbol{r})
\\

\end{align*}
$$

となります。

ここで、もし重なり積分がゼロ

$$
\int \chi_A(\boldsymbol{r})\chi_B(\boldsymbol{r})d\boldsymbol{r} = 0
$$

あるいはもう少し（今考えたい）重なり積分っぽく書いて、水素原子$A$の位置を原点、水素原子$B$の位置を$\boldsymbol{R}$と書いて

$$
\int \chi_s(\boldsymbol{r})\chi_s(\boldsymbol{r} - \boldsymbol{R}
)d\boldsymbol{r} = 0
$$

とすると、これは結局二つの波動関数の「重なり」がゼロで、$\chi_s(\boldsymbol{r})\chi_s(\boldsymbol{r} - \boldsymbol{R})=0$と考えると、

$$
\begin{align*}
\rho(\boldsymbol{r}) 

&=
 \frac{1}{\sqrt{2(1 + S)}}
\left(\chi_A^2(\boldsymbol{r}) + \chi_B^2(\boldsymbol{r})
\right)

+
\sqrt{\frac{2}{1 + S}}
 \chi_A(\boldsymbol{r}) \chi_B(\boldsymbol{r})
\\

&=
 \frac{1}{\sqrt{2}}
\left(\chi_s^2(\boldsymbol{r}) + \chi_s^2(\boldsymbol{r} - \boldsymbol{R})
\right)

\end{align*}
$$

となります。これは、ベクトル$\boldsymbol{R}$だけ離れた完全に孤立している水素原子（の電子）の分布を表しており、例えば距離が十分遠くに離れていた場合に対応していそうです。

一方$S\neq0$では、

$$
\begin{align*}
\rho(\boldsymbol{r}) 

&=
 \frac{1}{\sqrt{2(1 + S)}}
\left(\chi_A^2(\boldsymbol{r}) + \chi_B^2(\boldsymbol{r})
\right)

+
\sqrt{\frac{2}{1 + S}}
 \chi_A(\boldsymbol{r}) \chi_B(\boldsymbol{r})
\\


&=
 \frac{1}{\sqrt{2(1 + S)}}
\left(\chi_s^2(\boldsymbol{r}) + \chi_s^2(\boldsymbol{r}-\boldsymbol{R})
\right)

+
\sqrt{\frac{2}{1 + S}}
 \chi_s(\boldsymbol{r}) \chi_s(\boldsymbol{r} - \boldsymbol{R})

\\



\end{align*}
$$

です。これは、「重なり積分がゼロだった＝距離的に離れていた状態から、重なり積分が値を持つくらいまで近づいた際に、それぞれの水素原子が孤立していた場合の分布から

$$
\frac{1}{\sqrt{2}}- \frac{1}{\sqrt{2(1+S)}} 
$$

くらい確率分布が「減り」、代わりに$\chi_s(\boldsymbol{r}) \chi_s(\boldsymbol{r} - \boldsymbol{R})$の分の確率分布が$\sqrt{\frac{2}{1 + S}}$だけ「増えた」というように見ることができます。$\chi_s(\boldsymbol{r}) \chi_s(\boldsymbol{r} - \boldsymbol{R})$が大きな値を持つのは、両方の電子分布が小さくなりすぎていない分子の中間くらいなので、これは結局、2つの水素原子の間くらいに電荷分布が増え、その増えた電荷分布が原子間を繋ぎとめる「結合」の役割を果たすのだ、と解釈できます。

これはまさに結合性軌道と呼ぶにふさわしい軌道だと言えるわけです。

と、量子化学で分子の結合のあたりでなされる説明を雑に追ってきました。つまりここでは、二つの原子軌道間の重なり積分

$$
S = \int \chi_s(\boldsymbol{r})\chi_s(\boldsymbol{r} - \boldsymbol{R}
)d\boldsymbol{r} 
$$

が大きければ、結合性軌道

$$
\begin{align*}
\psi(\boldsymbol{r}) &= C(\chi_A(\boldsymbol{r}) + \chi_B(\boldsymbol{r}))\\

&=
\frac{1}{\sqrt{2(1 + S)}}
(\chi_A(\boldsymbol{r}) + \chi_B(\boldsymbol{r}))

\end{align*}
$$

は分子を形成し、小さければ形成しない（「中間地点」に電子が増えず、結合が形成されない）ものとして重なり積分を捉えることができそうです。（本当か？）

ただ、このような考え方は量子化学のように分子の結合を考える際に、「重なり積分が大きい（小さい）」＝「中間地点に電子のいる確率が大きく（小さく）、分子を作りやすい（にくい）」のように「重なり積分が大きければ/小さければどのようなことが起きるか」を考えるヒントにはなるのですが、「重なり積分」という量は何を表しているのか？」「なぜ重なり積分を小さいと置いてよいのか？」「重なり積分を小さいと置くことは何を意味しているのか？」について考えるには、まだ足りていないような気がします。そこで少し違う観点からこの積分を見ていこうと思います。

## 妄想

### 量子力学の公理（抜粋）

重なり積分の（ただ「波動関数の重なり」という以上の）物理的な意味について考える出発点として、これまでさらっと触れてきた量子力学の公理：

:::message

波動関数$\psi(\boldsymbol{r})$を、ある物理量$A$に対応する演算子$\hat{A}$の固有値$a_i$に対応する正規直交基底をなす固有関数$\phi_i(\boldsymbol{r}): \hat{A}\psi(\boldsymbol{r}) = a_i\phi_i(\boldsymbol{r})$で展開したとき：

$$
\psi(\boldsymbol{r}) = \sum_i c_i \phi_i(\boldsymbol{r})
$$

この波動関数$\psi(\boldsymbol{r})$に対して物理量$A$を測定した際、固有値$A_j$が得られる確率は

$$
|c_j|^2 = 
\left|\int \phi_j^*(\boldsymbol{r})\psi(\boldsymbol{r})d\boldsymbol{r}
\right|^2
$$

で得られる。

また、上記のように測定を行った後は、電子の波動関数はその測定で得られた固有値（物理量）に対応する固有関数で表される状態となる
:::

を考えることにします。

また、ここで

:::message
物理量$A$に対応する演算子$\hat{A}$の固有状態$\phi_i(\boldsymbol{r})$は正規直交基底をなす完全系である
:::

ことを仮定します。

この仮定はどの程度正しいのかは正直なところ私にはよくわかっていないのですが、J.J.サクライには「これは要請である」とか書いてあるので、まあとりあえず正しいことにします。（量子力学の公理系をどのように取るのかにもよるのかもしれませんが、あまり深いことは考えないことにして。。。）

### 重なり積分の物理的意味

というわけで、今考えたい重なり積分

$$
s_{\boldsymbol{0},\boldsymbol{R}}^{m,l} = \int \phi_m^*(\boldsymbol{r})\phi_l(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}
$$

に対応させて考えると、原点に中心を持つ原子軌道が固有状態となるような物理量があるような空間設定を考えれば、「格子点$\boldsymbol{R}$中心を持つ原子軌道と原点に中心を持つ原子軌道の（複素共役の）積の積分」に、「原点に中心を持つ原子軌道に対応した固有値が測定される確率」という意味を持たせられそうです。
最も手っ取り早い物理量はエネルギーだろうと思います。つまり、原点に、孤立した原子や分子ポテンシャルが一つだけあるような設定を考えれば、その系で測定されたエネルギーが、その孤立原子のハミルトニアンの固有状態と対応します。

さて、というわけで以下のような状況を考えてみます。

- 空間内には原点に孤立した原子核（又は内殻電子の平均場も含んだイオン）が存在し、この系のハミルトニアンを孤立原子核のポテンシャルを$V(\boldsymbol{r})$として$\hat{H} = (-\hbar^2/2m) \nabla^2 + V(\boldsymbol{r})$、このハミルトニアンの固有値$\varepsilon_m$を持つ固有関数を$\phi_m(\boldsymbol{r}): \hat{H}\phi_m(\boldsymbol{r}) = \varepsilon_m\phi_m(\boldsymbol{r})$とする。
- そのような空間内に、波動関数$\phi_l(\boldsymbol{r} - \boldsymbol{R})$で表される電子、つまり原点からベクトル$\boldsymbol{R}$離れた位置で$l$番目の原子準位を持つ原子軌道の電子が、（理由は分からないが）ポツンと浮かんでいる。

この時、仮定よりポツンと浮かんだ電子の波動関数$\phi_l(\boldsymbol{r} - \boldsymbol{R}$は、空間のハミルトニアンの固有状態$\phi_m(\boldsymbol{r})$を用いて

$$
\phi_l(\boldsymbol{r}-\boldsymbol{R}) = \sum_m c_m  \phi_m(\boldsymbol{r})
$$

と展開できるはずです。（本当なのだろうか？）

そしてこの時の展開係数$c_m$は、

$$
c_m = \int \phi_m^*(\boldsymbol{r})\phi_l(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r} = s_{\boldsymbol{0},\boldsymbol{R}}^{m,l}
$$

と重なり積分で表されることとなります。

量子力学の公理に従って考えると、ある波動関数を物理量の演算子＝ハミルトニアンの固有関数で展開した際の展開係数（の絶対値の2乗）は、その波動関数を観測した際、係数に対応する固有値が得られ、電子の波動関数がその固有値に対応する固有関数に「収束する」確率ですので、従って「重なり積分$s_{\boldsymbol{0},\boldsymbol{R}}^{m,l}$は

「原点に孤立原子ポテンシャルが存在する空間で、原点から$\boldsymbol{R}$離れた位置に（何故か）波動関数$\phi_l(\boldsymbol{r}-\boldsymbol{R})$が存在しているような状況において、電子のエネルギーを測定したら「原点を中心とする、準位$m$の原子軌道の固有値が得られる確率」

あるいはもう少し勝手な解釈のもと言い換えれば、

**「原点に孤立原子ポテンシャルが存在する空間で、原点から$\boldsymbol{R}$離れた位置に（何故か）波動関数$\phi_l(\boldsymbol{r}-\boldsymbol{R})$が存在しているような状況において、電子の状態を測定してみたら原点を中心とする原子準位$m$の原子軌道$\phi_m(\boldsymbol{r})$になっていた確率」**

を意味する積分であると言えるのではないでしょうか。

特に前章で考えた飛び移り積分との違いは、飛び移り積分は「原点に孤立した原子準位の**微小時間後**の位置（状態）」が隣の格子点に「飛び移っていた」確率を考えていました。飛び移り積分の形

$$
\int \phi_{m'}(\boldsymbol{r}-\boldsymbol{R}')V(\boldsymbol{r} - \boldsymbol{R}')\phi_m(\boldsymbol{r})d\boldsymbol{r}
$$

からもなんとなくイメージできるように、「格子点$\boldsymbol{R}$のポテンシャルに原点の電子が引っ張られているような感じです。

一方重なり積分の場合は、確かにハミルトニアンとしては孤立原子のハミルトニアンを考えており、孤立原子のポテンシャルは存在していますが、測定の時点では原点から$\boldsymbol{R}$離れた位置にいる準位$l$の電子$\phi_l(\boldsymbol{r} - \boldsymbol{R})$は、ポテンシャルを感じる暇もなく、原点の準位$m$の原子軌道$\phi_m(\boldsymbol{r})$になってしまっています。

## 重なり積分が小さい理由（の妄想）

以上で重なり積分のイメージがついたところで、最後に重なり積分を小さくしてよい（特に飛び移り積分は残して、重なり積分だけゼロとしてよい）事について納得するための妄想を繰り広げていきたいと思います。

ここから記号をそろえて、原点に局在した準位$m$の原子軌道と、格子点$\boldsymbol{R}$に局在した準位$l$の原子軌道$\phi_l(\boldsymbol{r} - \boldsymbol{R})$を考え、それらの間の飛び移り積分

$$
t_{(m,\boldsymbol{0})\leftarrow (l,\boldsymbol{R})} = 

\int \phi_m^*(\boldsymbol{r} )V(\boldsymbol{r}-\boldsymbol{R})\phi_l(\boldsymbol{r}-\boldsymbol{R})d\boldsymbol{r}
$$

と重なり積分

$$
s_{\boldsymbol{0},\boldsymbol{R}}^{m,l} = \int \phi_m^*(\boldsymbol{r})\phi_l(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}
$$

のそれぞれが持つ物理的 ~~妄想~~ イメージを比べてみましょう。どちらもTight-bindingモデルを考える際に出てくる積分ですが、前者は計算過程で残され、後者は無視されます。

前者については、前章で述べたように周期的な孤立原子のポテンシャルがあるような空間内で、格子点$\boldsymbol{R}$にいる原子準位$\phi_l(\boldsymbol{r} - \boldsymbol{R})$が原点の格子点に$\phi_m(\boldsymbol{r})$で飛び移る確率でした。

一方、後者の重なり積分は原点に一つ孤立原子のポテンシャルがあるような状況で、格子点$\boldsymbol{R}$にいる電子$\phi_l(\boldsymbol{r} - \boldsymbol{R})$を測定してみたら、原点の孤立原子の電子軌道$\phi_m(\boldsymbol{r})$になってる（？）確率でした。

前者は格子点のポテンシャルに束縛された状態から別の格子点に飛び移る確率であるのに対し、後者は何にも縛られないいわばふわふわ浮かんだ原子軌道が、原点の孤立原子に束縛される確率を表しており（と妄想しており）、一見後者の方がありえそうに思える気もします。

しかしよく考えてみると、格子点$\boldsymbol{R}$に「ふわふわ漂う」原子軌道と原点の孤立原子ポテンシャルがあるとしましょう。原子軌道というのも結局は、中心力の引力ポテンシャルに束縛された電子状態なので、まあ調和振動子、あるいは天体の周りをまわる惑星のような運動をしているのでしょう（実際のところ原子軌道の電子がどのような運動をしているのか知るすべは（今のところ）無いのだと思いますが）

というわけで電子がブランコを漕いでいる少年のようなイメージをしておきましょう。

![](/images/tb/overlap-1.png)

さらに、格子点$\boldsymbol{R}$にいる「原子軌道」は特にポテンシャルに束縛されているわけでもなく「ふわふわ漂って」いるわけなので、上記のようなブランコ少年の、ブランコが突然消滅したような感じだと思われます。（ほんまに？）

![](/images/tb/overlap-2.png)


そしたら、ブランコ少年（原子軌道）はどこか遠くに飛んでいくこともあり得るし、

![](/images/tb/overlap-3.png)

もとの準位より高かったり、低かったりするエネルギー順位で孤立原子に束縛されることもあり得そうです。

![](/images/tb/overlap-4.png)

つまり様々な可能性があり得る分、その可能性のうちの一つ「格子点$\boldsymbol{R}$にいる電子$\phi_l(\boldsymbol{r} - \boldsymbol{R})$を測定してみたら、原点の孤立原子の電子軌道$\phi_m(\boldsymbol{r})$になってる（？）確率」は小さいことが納得できそうです。

もうちょい数学的（？）に考えると、原点からベクトル$\boldsymbol{R}$平行移動された原子軌道$\phi_l(\boldsymbol{r} - \boldsymbol{R})$を原点を中心とする原子軌道$\{\phi_m(\boldsymbol{r})\}$で展開したとき、前者と後者は中心がずれたあまり似ていない関数なので、後者のうちどれか一つの係数が大きくなるようなこともなく、均一な加重で展開されることになりそうです。

一方飛び移り積分の場合は、微小時間後の波動関数つまり原点の原子軌道へハミルトニアンを作用させた関数の展開

$$
\begin{align*}
\psi(\boldsymbol{r},\Delta t) &\simeq 

\psi(\boldsymbol{r},0 ) + \frac{-i}{\hbar}\Delta t\hat{H}\psi(\boldsymbol{r},0)
\\

&\simeq
\phi_m(\boldsymbol{r})
+
\frac{-i}{\hbar}\Delta t\hat{H}\phi_m(\boldsymbol{r})
\\

&\simeq

\left\{1-\frac{-i}{\hbar}\Delta t
(\varepsilon_m^{\rm a} + \Delta\varepsilon_{mm})
\right\}
\phi_m(\boldsymbol{r})\\

&\>\>\>\>+
\frac{-i}{\hbar}\Delta t\sum_{m'} \Delta\varepsilon_{m'm}\phi_{m'}(\boldsymbol{r})\\

&\>\>\>\>+\frac{-i}{\hbar}\Delta t\sum_{\boldsymbol{R}'}\sum_{m'}\left(-t_{\boldsymbol{R}'}^{mm'}\right)^*\phi_{m'}(\boldsymbol{r} + \boldsymbol{R}')


\end{align*}
$$

が示すように、原子軌道の「行先」は各格子点しかなく、特に最も近い最近接格子への飛び移りの確率は、（どこに行くかわからない重なり積分＝ふわふわ漂った原子軌道がある原子軌道に一致する確率に比べれば）それなりに大きくなりそうであろうと思われるわけです。





# おわりに

好き勝手なことを書いた本章にまとめも何もないような気がしますが、一応得られた結論として、異なる原子軌道・格子点上の重なり積分


$$
s_{\boldsymbol{0},\boldsymbol{R}}^{m,l} = \int \phi_m^*(\boldsymbol{r})\phi_l(\boldsymbol{r} - \boldsymbol{R})d\boldsymbol{r}
$$

は、「原点に孤立原子のポテンシャルがある状況において、格子点$\boldsymbol{R}$にたまたま原子軌道$\phi_l(\boldsymbol{r} - \boldsymbol{R})$がいた時、その電子を測定してみると原点の孤立原子の固有状態$\phi_m(\boldsymbol{r})$になった確率」のようなもんだということが分かりました。なんじゃそら。

そしてよくわかりませんが、この確率は「周期的に並んだ局所ポテンシャルの内どれか、特に再隣接の格子点のポテンシャルに飛び移る確率」よりは十分小さそうですね。というわけで我々は安心して、これからも異なる軌道・格子点間の重なり積分を（飛び移り積分と比べて）ゼロと近似していきましょう。


