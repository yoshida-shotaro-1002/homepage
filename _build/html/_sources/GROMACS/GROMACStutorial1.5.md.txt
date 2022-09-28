# 5.【補足】mdpファイルの詳細

ここでは、タンパク質のシミュレーションで使用したmdpファイルの詳細に触れる。 MDシミュレーションの一般的な流れと使用したmdpファイルは（ions.mdpを除いて）次のようなものであった。

1. エネルギー最小化 : minim.mdp
2. NVT平衡化 : nvt.mdp
3. NPT平衡化 : npt.mdp
4. 本計算 : md.mdp

それぞれのファイルでどのような設定を行っているのか、重要なポイントに絞って解説する。 ここで紹介しなかったコマンドやより詳しい解説は[GROMACS.mdpオプション](https://shalad2.github.io/JellyPage/gromacs/gromacs_mdp.html)または[公式ドキュメント](https://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html)を参照してほしい。 なお、全ての設定に対してコメントが付いているので、それを読んで理解してもらってもよいだろう。

## minim.mdp

エネルギー最小化のために使用したminim.mdpの内容は次のようなものであった。

```
minim.mdp
; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform
    
; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
```

1つ目のセクションでは、エネルギー最小化のために`integrator = steep`で最急降下アルゴリズムを採用している。 `emtol = 1000.0`は最大力Fmaxが1000 kJ mol-1 nm-1以下の値となるまでエネルギー最小化を行う設定をしている。 エネルギー最小化の最大ステップ数は`nsteps = 50000`としていて、これより前にエネルギーが最小化された場合はそこで停止する。

また、2つ目のセクションでは非結合相互作用（LJ相互作用、静電相互作用）に関する設定をしている。 `pbc = xyz`は3次元の周期的境界条件を指定している。

## nvt.mdp

300 Kにおける100 psのNVT平衡化のために使用したnvt.mdpの内容は次のようなものであった。

```
nvt.mdp
title                   = OPLS Lysozyme NVT equilibration 
define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
```

NVT平衡化のMD計算を行うために`integrator = md`としている。 タイムステップは2 fsに設定するので、psの単位で`dt = 0.002`としており、100 psのシミュレーションのためには`nsteps = 50000`が必要になる。

シミュレーションを始める際には原子に初速度が必要であるので、シミュレーション開始時に`gen_vel = yes`でマクスウェル分布に従ってランダムな速度を与える。 温度制御は`tcoupl = V-rescale`という方法で行い、温度は`ref_t = 300`を参照する。 NVTアンサンブルなので圧力は制御せず、`pcoupl = no`としている。

結果の出力は「Output control」の部分で行い、原子の座標や速度、系のエネルギー、シミュレーションのログファイルを何ステップおきに出力するかを設定する。 タイムステップは2 fsなので、`nstxout = 500`とすると座標の出力は1 psおきに設定される。

## npt.mdp

300 K、1barにおける100 psのNPT平衡化のために使用したnpt.mdpの内容は次のようなものであった。

```
npt.mdp
title                   = OPLS Lysozyme NPT equilibration 
define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes       ; Restarting after NVT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no        ; Velocity generation is off
```

NPT平衡化のMDの設定はNVTとほとんど同じであるが、NPTアンサンブルでは圧力制御を行うため`pcoupl = Parrinello-Rahman`としている。 この時の圧力の値は`ref_p = 1.0`を参照している。 温度制御はNVT平衡化の設定のままで継続して行っている。

またNPT平衡化では、系の原子に初速度を与える代わりにNVT平衡化の最終状態を引き継ぐため、`gen_vel = no`及び`continuation = yes`へ変更をしている。 こうすることで、指定したチェックポイントファイルを読み込み、シミュレーションを継続することが可能になる。

## md.mdp

最後に、本計算で用いたmd.mdpの内容を確認する。 本計算はNPTアンサンブルで行っているので、シミュレーション条件はNPT平衡化とほとんど同じである。

```
md.mdp
title                   = OPLS Lysozyme NPT equilibration 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 500000    ; 2 * 500000 = 1000 ps (1 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System    ; save the whole system
; Bond parameters
continuation            = yes       ; Restarting after NPT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel                 = no        ; Velocity generation is off
```

1 nsのシミュレーションを行っているので`nsteps = 500000`が必要である。

本計算での設定で特徴的なのは、今まで指定していた`define = -DPOSRES`を指定していないことである。 これはタンパク質の位置拘束を有効にする設定で、平衡化の段階ではタンパク質が大きく動かないようにするために指定していた。 しかし、本計算ではタンパク質のダイナミクスを観測するために指定を外している。

また出力形式であるが、今回は座標のみを書き出したいので、ファイルサイズの小さいxtc形式のトラジェクトリを指定している。 そのためtrr形式は`nstxout = 0`として出力せず、`nstxout-compressed = 5000`として5000ステップおきにxtc形式のトラジェクトリを出力している。 これで10 psおきに座標が記録されることになるので、1 nsのシミュレーションでは約100ステップ分のトラジェクトリを得ることになる。 VMDでトラジェクトリを眺める時に、読み込まれたステップ数を確認してほしい。

## まとめ

mdpファイルの形式と、その大まかな内容を見てきた。 詳細な説明はしていないので、次のリンクを参考にしてさらに理解を深めてほしい。

1. [GROMACS.mdpオプション](https://shalad2.github.io/JellyPage/gromacs/gromacs_mdp.html)
2. [公式ドキュメント（外部ページ）](https://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html)
