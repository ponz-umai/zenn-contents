---
title: "z.comがサービス終了するのでGWSに移行しようと思ったら、、、結局お名前メールと無料アカウントのGoogleサイトで足りた話（中編）"
emoji: "✈️"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [google,googlecloud,web]
published: true
---
## はじめに　前回の振り返り

[前編](https://zenn.dev/ponzumai/articles/blog-migration-memo)では、独自ドメインを使った最低コストのWEBページ公開と、メールの立ち上げ方法を模索した結果、当初はGWSを独自ドメインで契約したうえでGoogleサイト・Gmailを使う方針で進めていたのですが、色々設定していくうちに無料のGoogleアカウントで作成したGoogleサイトに独自ドメインを設定できることに気づいてしまい、「無料のGoogleアカウント＋ドメイン管理サービスが提供するメールサーバ」の運用に行きつきました。

てなわけでドメインを管理するお名前ドットコムが提供するサービスお名前メールにてメールサーバを立ち上げようとしていたわけですが、、、、謎のエラーでお名前メールに申し込めまず、時間切れで一旦終了したわけでした。

あれから3日くらいたっても一向にエラーが治らないし、お名前ドットコムに問い合わせを送っても返ってこないし、もうどないやねんということで違う選択肢を探すことにしました。

## 新たな方針：スタードメイン

で、もう別のドメイン管理サービスのメールサーバを契約するかと色々探していたところ、なんとドメイン登録費だけでレンタルサーバも、メールサーバも、サイト転送も使えるというサービス[スタードメイン](https://www.star-domain.jp/)を見かけまして、まあそれなりに歴史のある怪しくないところらしいと思えたのでこちらに全部乗り換え、メールもこっちで運用することにしました。

### まずは1ドメインを移管してテスト

とはいえ一気にやるのも怖いし、しかしZ.comのサービス終了期限はもはや目の前に迫っているので、丁度Z.comで取得していたドメインをスタードメインに移管して諸々使い勝手を試してみることにします。

移管申請から、メールで承認して、もう一回承認して、、、なんやかんややってると2～3日くらいで移管完了しました。
そこから次は無料サーバ利用申し込みをして、それも1日くらい待って、ついにドメイン移管と無料メールが使えるようになりました！

![alt text](/images/image-3.png)

### Googleサイトに独自ドメインを紐づけ

[こちらのブログ](https://www.value-domain.com/media/domain-usage-homepage/)の内容ほぼその通りに実行していくことで、（多少の待ち時間、1時間以内くらい）で独自ドメインの設定が完了します。素晴らしい・・・！

### メールアカウントの発行

![alt text](/images/image-4.png)メールアカウント発行はサーバ管理画面「メールアドレス設定」から、必要事項入力して登録ボタンを押せばすぐに使えるようになりました。簡単すぎでは？
WEBメールツールから実際にメールも送受信できます。

やばくない？

もうこれで十分なので、本命のドメインも~~お問い合わせが返ってこないし自社のサーバも管理できない~~お名前ドットコムから移管していくことにします。なんてこった。

既にGWSに紐づけてしまってますが、とりあえず何も考えずに移管手続きをしてみます。DNSとかよくわかってないのですが多分何とかなるのでしょう。多分。。。GWSの無料期間があと1週間くらいなので、できればそれまでに移管を完了して解約したい！

## Googleサイト運用のあれこれ

~~さて、無事Googleサイトへの独自ドメインの紐づけは完了したのですが、元々運用していたサイトが「www」無しで公開していたからなのか（Googleサイトは強制的にwwwをつけられる）、また別の原因なのか、直接www付きでURLにアクセスしたら上手く表示されるのですが、www無しでアクセスしてみると、↓のようにエラーが出ちゃいます。~~

![alt text](/images/image-5.png)

~~証明書を見てみるとスタードメインのSSL証明書が適用されている？みたいなのですが、、、マジで知識がないので何もわかりません。~~

~~これでいいのかどうか何もわからないのですが、スタードメインはサイト転送設定も無料でできるのでwww無しを、ありに転送するような設定をしておきました。（これで合ってるのか？）~~

~~まあ、エラー文にもしばらく時間をおいて、と書いてあるので、ここは数日様子を見てみるか、暇なときにググってみるかしてみようと思います。根本的にはDNSとか、なんか色々勉強したらわかるのだろうけど。スタードメインはSSL設定も無料でできるので（私は回し者か？）、設定してみたら上手く行くのか～？
まあなんか数日おきに適当にやってみようと思います。~~

（2024年7月10日追記　全て適当なことを書いていました、詳細は後編へ）


## おわり

ひとまず、最大の目標だったZ.com サービス終了までに独自ドメインのサイトと、メールアカウントの移行方針の目途はついたのでめでたいところ。
あと最後、後編にGWSの解約と、謎のエラー解決について書いてみようと思います。
