# ArtemHopper 
# Calculations by Jarvist Moore Frost, 2018 
# These codes were developed with Julia 0.5.0, and requires the Optim and Plots packages.

# This file, when run under Julia, should regenerate polaron data associated with Tom Hopper / Artem Bakulin's paper on cooling of Perovskites. 

push!(LOAD_PATH,"../PolaronMobility.jl/src/") # load module from local directory

using PolaronMobility 
using PlotPolaron # Plots dependency

##### load in library routines... #####
# Plot figures with Plots, which defaults to Pyplot backend
#using Plots
#pyplot() # PyPlot (matplotlib) backend, to Plots
#gr() # GR backend, to Plots
#default(grid=false) # No silly dotted grid lines
#default(size=(400,300)) # A good small size for two-column EPS output
#default(size=(800,600)) # Nice size for small-ish PNGs for slides

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space
# Units
Å=1E-10
const THzInCm1=0.02998

### Frequency + IR activity data; extracted from VASP via vasp_ir.sh

mutable struct PolarMaterial
    name
    IR
    ϵoptic; ϵstatic
    effectivemass
    effectivefrequency
end
PolarMaterial()=("",[],0,0,0,0) # default

# Cross-material defaults for consistency
IodideMass=0.12
BromideMass=1.2 * IodideMass

# ϵoptic - expt results from Leguy et al. Expt + Theory optical props
# 10.1039/C5NR05435D
# Supp. Info (.txt file) - first line 
#   'n' refract index for single crystals @ 1.22 eV
Iodideϵoptic=2.369146846^2 # MAPI ϵoptic = 5.612856777911748
Bromideϵoptic=2.041097522^2 # MAPBr ϵoptic = 4.16607909431454


### 2018-05
# /work/jmf02/2018-05-TomHopper-MAPIs/0002-InitialAttemptDielectric
#   CsPbBr.vasp  FAPbI.vasp  MAPbBr.vasp  MAPbI.vasp

# MAPbI/.../exact.res.txt
MAPbI=PolarMaterial("MAPbI",
[
001 3208.947236 .51670227836603600004
002 3106.662725 1.7640215248717381
003 3098.051235 .61587964854276
004 3084.528751 .0076969866838925
005 3077.667580 .0960736882203565
006 2981.869101 .005014457498120529
007 1563.922801 .25712100479294052165
008 1548.256737 .14676293861822640000
009 1468.853932 .216435672793777129
010 1430.780855 .013590026716589604
011 1423.354182 .029615732196994192
012 1383.244893 .017257965681922324
013 1237.053096 .001051164702816538
014 1218.542713 .0216564976644344
015 1007.405726 .008956671839704644
016 912.543801 .046965342791386625
017 878.882843 .055571480292890729
018 317.093051 .000743955374021561
019 133.227525 .085542805843646177
020 129.839102 .007671109865113476
021 117.379891 .058323329362874064
022 91.577572 .026910992923781225
023 81.014431 .234367016325789800
024 74.966229 .273155780653146369
025 68.974799 .215920083638021636
026 65.778507 .067005114863164404
027 58.353707 .0495305352455325
028 34.549682 .010939016415544804
029 32.735898 .00516853627954483249
030 32.596119 .009046068689192336
031 28.748393 .0136761595930616
032 26.680367 .00134674307868884324
033 18.668885 .00717601082780908416
#034 0.035663 .00000000571618602225
#035 0.049863 .00000003437152336349
#036 0.061231 .00000000835418188889
],
Iodideϵoptic, Iodideϵoptic+15.4, IodideMass, 0.0)

#> cat FAPbI/intensities/results/exact.res.txt 
FAPbI=PolarMaterial("FAPbI",
[
001 3471.688121 .27863219234217113609
002 3469.582090 1.2822244440277064
003 3379.633172 .0584829230948114
004 3342.478686 .52095193439234023580
005 3092.671091 .064941803812768756
006 1737.055377 .39089879800791732502
007 1596.565531 .3264853728772025
008 1517.171587 .00014256441352220823
009 1354.933816 .08645731489637966564
010 1322.477933 .27244760209040086954
011 1133.153739 .002545568132577949
012 997.855314 .01621445690730287801
013 989.981489 .02646413168677242851
014 678.959967 .00581198869711315460
015 572.933586 .00000005500153921618
016 511.491098 .91844239563011105816
017 503.832211 .01510179649495881876
018 456.510291 .00000380771407857007
019 90.705715 .025128707642230149
020 81.638692 .16374810390589997416
021 76.008939 .04675282657734450000
022 65.967822 .0119352927940765
023 65.156075 .362801902096225289
024 64.703260 .0247183893078562
025 45.359118 .27752364960532730749
026 36.838820 .01249030092244275488
027 32.888165 .00012814109833170609
028 29.641620 .02150653349713943234
029 26.147604 .004563107311737530
030 24.619597 .00418726688573374436
#031 0.389193 .00000177033868247888 # Acoustic
#032 0.580609 .00010047866332146178
#033 0.813049 .00000255093586843149
#034 8.391534 .08518006103778671742 # Soft modes - -ve ?
#035 48.796868 .02787431057622411400
#036 104.953267 .00000005953029888586
],
Iodideϵoptic, Iodideϵoptic+22.9, IodideMass, 0.0) 
# Not sure why ϵionic is larger than MA
# See jmf02@login-2-internal:/work/jmf02/2018-05-TomHopper-MAPIs/0002-InitialAttemptDielectric/ 
# Possibly poorly converged structure? --> overly soft modes

#> cat CsPbBr/intensities/results/exact.res.txt 
CsPbBr=PolarMaterial("CsPbBr",
[
001 102.584848 .298397802564
002 102.574769 .298414190529
003 102.546547 .298578094929
004 38.616236 .011199565584
005 38.604913 .011131938064
006 38.600490 .011180736121
007 19.008486 .00000004012009
008 19.007460 .000000037315421584
009 19.006706 .000000035940955561
#010 3.700461 .00088522315729 # acoustic
#011 3.697970 .00088594927201
#012 3.667602 .00087957916929
#013 20.032115 .023329813081 #soft
#014 20.095315 .023257775025
#015 20.098301 .023257775025
],
Bromideϵoptic,Bromideϵoptic+10.0, BromideMass,0.0)


#> cat MAPbBr/intensities/results/exact.res.txt 
MAPbBr=PolarMaterial("MAPbBr",
[
001 3162.380018 1.26914271741417
002 3148.741031 .44683916306084
003 3102.946830 .885064155116179409
004 3084.321691 .064698458045588769
005 3080.836830 .00095582401515090400
006 2983.348944 .009364378765001961
007 1577.634841 .22505651278592845409
008 1545.090612 .15517853860628
009 1489.174432 .177091072101164736
010 1431.355566 .029118040066118961
011 1423.999303 .024060432348775746
012 1388.647856 .016612483524033801
013 1243.896741 .00004389292567414004
014 1222.263661 .00532165455936914249
015 1017.122520 .00625838293840781561
016 918.684312 .036019157350986260
017 890.816087 .022626892583943696
018 328.676462 .00908628023356039930
019 159.092720 .03411041882947896096
020 150.197776 .078090686403216329
021 109.570391 .0118428813453194
022 105.920699 .0834580502620736
023 92.909502 .16878545309356041041
024 83.044727 .321103336416702164
025 74.292980 .03026873239459394384
026 72.897653 .3597632983323201
027 62.534207 .0270758068177441
028 47.592144 .00434340091536630569
029 42.857313 .006538777905991889
030 41.547462 .0007395162249734
031 39.716407 .01750481976953124811
032 29.156555 .001542205890170336
033 18.902352 .01209986847211389125
#034 0.185557 .00000000972650667921 # acoustic
#035 0.239986 .00000000642548283024
#036 0.419451 .00000026687618127000
],
Bromideϵoptic,Bromideϵoptic+15.9, BromideMass,0.0)

# Re-relaxed structure
# CX1: /work/jmf02/2018-05-TomHopper-MAPIs/0004-DielectFromRelaxed
FAPbBr=PolarMaterial("FAPbBr",[
001 3495.800077 1.3546770943165760
002 3488.457140 .1492227939197984
003 3393.771847 .134401358853954821
004 3353.723320 1.1285076514103970
005 3101.150462 .02263645613873268569
006 1745.176076 .918746227323340461
007 1604.109503 .265035531778307325
008 1506.342749 .047894393119858901
009 1358.792791 .01564042871697359621
010 1321.446232 .10305105674046726244
011 1158.559509 .0007527932918342
012 994.948324 .0772714383483509
013 985.926164 .06621811883430
014 680.202054 .0225777162104969
015 576.154307 .0003346695850566
016 508.434617 .001936195824940889
017 491.507527 .9130197360644597
018 398.048840 .0001304997138437
019 114.125230 .03978733406892744016
020 100.790513 .223458531117122481
021 89.353095 .03942932209685728144
022 88.341500 .02225158276393
023 68.995007 .411348594568822525
024 50.642030 .016380246970744549
025 47.598506 .0035914665924504
026 45.975085 .1475779565907290
027 41.604417 .01714947753012338861
028 37.500633 .0093255598013384
029 35.047638 .0116313388610132
030 33.722222 .006948027332700041
031 18.480578 .1793392652005945
032 0.891908 .000001229792292236
033 1.112119 .000001042480412525
#034 1.550482 .000016948559987016 # acoustic
#035 36.640240 .018664093936020724
#036 175.550538 .0000641152821203
],
Bromideϵoptic,Bromideϵoptic+33.79612115611553, BromideMass,0.0)


"""
E=25E-3*q # 25 meV, in Joules
println("Test for k_B T=E J; i.e. thermal energy phonon occupation at STP : ",BE(300,E))
"""
BE(T,En)=1./(exp.(En/(kB*T)) - 1.0) # Bose-Einstein Occupation

function phononinfo(PhononFreqs; T=300)
    PhononEnergy = PhononFreqs.*1E12* 2*π * ħ # simply ħω ; including unit conversions
    println("Phonon energies, in meV: ",1000*PhononEnergy/q)
    println("Sum of phonon energies, in meV: ",sum(1000*PhononEnergy/q))
    
    Z=sum(exp.(-PhononEnergy/(kB*T)))
    # Now, we've just calculated the partition function (Z) for a single unit cell, at a certain temperature

    # Boltzmann / Thermodynamic occupation
    #ThermalOccupation=exp.(-PhononEnergy/(kB*T))/Z
    
    println("BE thermal occupation of phonon states...")
    BEThermalOccupation=BE(T,PhononEnergy)
    println(BEThermalOccupation)

    #println("(Average Thermal) Energy in phonons (BE stats) of unitcell at T= $T, ", 1000/q*BEThermalOccupation.*PhononEnergy, " meV")
    println("Sum: ",sum(1000/q*BEThermalOccupation.*PhononEnergy)," meV")
end

for material in [FAPbI, FAPbBr,  MAPbBr, MAPbI, CsPbBr ]
    println("\n\nMaterial: ",material.name)

    freq=material.IR[:,2]*THzInCm1
    activity=material.IR[:,3]
    
    hellwarth=[freq activity]
    if length(hellwarth)>18
        hellwarth=hellwarth[18:end,:] # if hybrid organic; delete intra-molecular modes
        println(hellwarth)
    end
#    show(HellwarthBScheme(hellwarth))
    material.effectivefrequency=HellwarthBScheme(hellwarth)

    println("Effective frequency: ",material.effectivefrequency)

    mob=polaronmobility([300,1000], material.ϵoptic, material.ϵstatic, material.effectivefrequency*1E12, material.effectivemass)

	phononinfo(freq,T=300)
    println("Material: ",material.name, " α: ",mob.α, 
            " Scatterτ: ",mob.Tau," rfsi: ",mob.rfsi)
    mob.rfsi
    mob.α
end


end

# --> 0.12 for electrons and 0.15 for holes, in MAPI. See 2014 PRB.
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=
MAPIe=polaronmobility(10:20:1000, 4.5, 24.1, 2.25E12, 0.12)
savepolaron("MAPI-electron",MAPIe)
plotpolaron("MAPI-electron", MAPIe)
MAPIh=polaronmobility(10:20:1000, 4.5, 24.1, 2.25E12, 0.15)
plotpolaron("MAPI-hole", MAPIh)
savepolaron("MAPI-hole",MAPIh)

## Polaron radius vs. Temperature
# More complex figure for the Thermal Pathways paper;
# https://github.com/WMD-group/2017-03-ThermalPathways/
p=MAPIe 

plot(p.T,p.rfsi./Å, markersize=2,marker=:rect,
    label="Polaron radius",xlab="Temperature (K)",ylab="Polaron Radius (Angstrom)",ylims=(0,Inf))
plot!(p.T,p.rfsmallalpha./Å,label="T=0 Schultz small alpha polaron radius")
savefig("ArtemHopper-radius.png")
savefig("ArtemHopper-radius.pdf")

plot(p.T,p.v./p.w,label="v/w",markersize=2,marker=:circle,
    xlab="Temperature (K)", ylab="\hbar\omega")
savefig("ArtemHopper-vwratio.png")
savefig("ThermalPathways-vwratio.pdf")

    plot!(p.T,p.v,label="v",markersize=2,marker=:uptriangle)
plot!(p.T,p.w,label="w",markersize=2,marker=:downtriangle)
savefig("ArtemHopper-vw.png")
savefig("ArtemHopper-vw.pdf")

println("That's me!")

