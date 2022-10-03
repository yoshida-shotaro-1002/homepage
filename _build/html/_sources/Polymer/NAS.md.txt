# NASの使い方

## NASって何？

デスクトップPCだけだと大規模計算のデータを保存するのには足りないため、大容量の記憶媒体であるNASをネットワーク中でマウント接続して保存することで大規模なデータを管理している。

## NASをマウントする

マウントするディレクトリを以下のように作成する。

```shell
sudo mkdir /mnt/ImPACT01
sudo mkdir /mnt/ImPACT02
sudo mkdir /mnt/point01
sudo mkdir /mnt/point02
```
次に以下のコマンドを打ち込み、以下の画面のように下 5 行を追記する。必要に応じて頭にsudoをつける。

```shell
vi /etc/fstab
```

```eval_rst
.. figure:: image/fstab.png
    :width: 80%
```

編集できたら保存して閉じ、以下のコマンドを打ち込んでマウント接続を実行する。

```shell
sudo mount -a
```
実行後、ImPACT02のディレクトリに入り、`ls`を実行してマウントされているか確認する。
