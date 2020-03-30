# by Héctor Otárola and Rodrigo Moreno
using JuMP
using Cbc, CSV

#Change to folder path
cd("/home/script/")

m = model = JuMP.Model(with_optimizer(Cbc.Optimizer, logLevel=1))

VecPar = CSV.read("VecPar.csv", header = false)
Dem = CSV.read("Dem.csv", header =false)[:,1]
PVprofile = CSV.read("PVprofile.csv", header = false)[:,1]
Wprofile = CSV.read("Wprofile.csv", header =false)[:,1]
PreFacPV = CSV.read("PreFacPV.csv", header = false)
PreFacB = CSV.read("PreFacB.csv", header = false)


Rcost = VecPar[1,1]
CInvCon = VecPar[2,1]
CInvWind = VecPar[3,1]
CInvPV = VecPar[4,1]
CInvBat = VecPar[5,1]
CInvPS = VecPar[6,1]
CInvLin = VecPar[7,1]
COpeCon = VecPar[8,1]
COpeCoal = VecPar[9,1]
COpeGas = VecPar[10,1]
COpeOil = VecPar[11,1]
CEncCoal = VecPar[12,1]
CApaCoal = VecPar[13,1]
CoalCap = VecPar[14,1]
GasCap = VecPar[15,1]
OilCap = VecPar[16,1]
MinTecCoal = VecPar[17,1]
RampCoal = VecPar[18,1]
RampGas = VecPar[19,1]
RampOil = VecPar[20,1]
DurStor = VecPar[21,1]
EficiencyStor = VecPar[22,1]
DemBase = VecPar[23,1]
RDem = VecPar[24,1]
LinCap_12 = VecPar[25,1]
LinCap_23 = VecPar[26,1]
LinCap_23n = VecPar[27,1]
LinCap_13 = VecPar[28,1]
XLin = VecPar[29,1]
PSCap = VecPar[30,1]

nY = 3 #Amount of years
nS = 4 #Amount of nodes of scenario-tree nodes
nT = 24 #Amount of hours
BigM = 10000

@variables m begin
0 <= Pc[1:nT,1:nY,1:nS] <= CoalCap
0 <= Pg[1:nT,1:nY,1:nS] <= GasCap
0 <= Po[1:nT,1:nY,1:nS] <= OilCap
Xc1[1:nT,1:nY,1:nS], Bin
Ec1[1:nT,1:nY,1:nS], Bin
Ac1[1:nT,1:nY,1:nS], Bin
Pc2[1:nT,1:nY,1:nS] >= 0
W[1:nT,1:nY,1:nS] >= 0
PV[1:nT,1:nY,1:nS] >= 0
BW[1:nT,1:nY,1:nS]
BWc[1:nT,1:nY,1:nS] >= 0
BWd[1:nT,1:nY,1:nS] >= 0
BWe[1:nT,1:nY,1:nS] >= 0
Ps[1:nT,1:nY,1:nS]
BPVe[1:nT,1:nY,1:nS] >= 0
BPVd[1:nT,1:nY,1:nS] >= 0
BPVc[1:nT,1:nY,1:nS] >= 0
BPV[1:nT,1:nY,1:nS]

F12[1:nT,1:nY,1:nS]; F23[1:nT,1:nY,1:nS]; F23n[1:nT,1:nY,1:nS]; F13[1:nT,1:nY,1:nS]; T1[1:nT,1:nY,1:nS]
T2[1:nT,1:nY,1:nS]; T3[1:nT,1:nY,1:nS]

WMax[1:nY,1:nS] >= 0
PVMax[1:nY,1:nS] >= 0
BWMax[1:nY,1:nS] >= 0
BPVMax[1:nY,1:nS] >= 0
Pc2_max[1:nY,1:nS] >= 0
XPS[1:nY,1:nS],Bin
XFCo2[1:nY,1:nS],Bin
end

@expression(m, An_Inv[y=1:nY,s=1:nS], CInvCon*Pc2_max[y,s] + CInvWind*WMax[y,s] + PreFacPV[y,s]*CInvPV*PVMax[y,s] + PreFacB[y,s]*CInvBat*BWMax[y,s] + CInvLin*XFCo2[y,s] + BPVMax[y,s]*PreFacB[y,s]*CInvBat*10000 + CInvPS*XPS[y,s])
@expression(m, An_Op[y=1:nY,s=1:nS],  365*(sum(COpeCon*Pc2[t,y,s] + COpeCoal*Pc[t,y,s] + COpeGas*Pg[t,y,s] + COpeOil*Po[t,y,s] + CEncCoal*Ec1[t] + CApaCoal*Ac1[t]  for t=1:nT)))

@objective(m, Min, sum(0.25*sum((An_Inv[y,s]+An_Op[y,s])/(1+Rcost)^(5*(y-1)) for y=1:nY) for s = 1:nS))


for t = 1:nT
  for y = 1:nY
    for s = 1:nS
        @constraints(m,begin
        #Power balance constraints (Kirchhoffs current law)
          Pg[t,y,s] == F12[t,y,s] + F13[t,y,s]
          Pc[t,y,s] + Pc2[t,y,s] + W[t,y,s] + BW[t,y,s] + F12[t,y,s] == F23[t,y,s] + Ps[t,y,s] + F23n[t,y,s]
          Po[t,y,s] + PV[t,y,s] + BPV[t,y,s] + F13[t,y,s] + F23[t,y,s] + F23n[t,y,s] + Ps[t,y,s] == DemBase*Dem[t]*(1+RDem)^(5*(y-1))

    #Power transfers capacity constraints
          F12[t,y,s] <= LinCap_12
          F23[t,y,s] + Ps[t,y,s] <= LinCap_23
          F23n[t,y,s] <= LinCap_23n*XFCo2[y,s]
          F13[t,y,s] <= LinCap_13
          Ps[t,y,s] <= PSCap*XPS[y,s]


          F12[t,y,s] >= -LinCap_12
          F23[t,y,s] + Ps[t,y,s] >= -LinCap_23
          F23n[t,y,s] >= -LinCap_23n*XFCo2[y,s]
          F13[t,y,s] >= -LinCap_13
          Ps[t,y,s] >= -PSCap*XPS[y,s]

    #Power flow constraints (Kirchhoffs voltage law)
          F12[t,y,s] == (T1[t,y,s] - T2[t,y,s])/XLin
          F23[t,y,s] == (T2[t,y,s] - T3[t,y,s])/XLin
          F13[t,y,s] == (T1[t,y,s] - T3[t,y,s])/XLin
    #Big-M constraints for the new line
          F23n[t,y,s] <= (T2[t,y,s] - T3[t,y,s])/XLin + BigM*(1-XFCo2[y,s])
          F23n[t,y,s] >= (T2[t,y,s] - T3[t,y,s])/XLin - BigM*(1-XFCo2[y,s])

          T1[t,y,s] == 0

    #Power output capacity constraints
          Pc[t,y,s] <= CoalCap*Xc1[t,y,s]
          Pc2[t,y,s] <= Pc2_max[y,s]
          Pc[t,y,s] >= MinTecCoal*CoalCap*Xc1[t,y,s]
          W[t,y,s] <= WMax[y,s]*Wprofile[t]
          PV[t,y,s] <= PVMax[y,s]*PVprofile[t]

    #Battery constraints
          BW[t,y,s] <= BWMax[y,s]
          BPV[t,y,s] <= BPVMax[y,s]
          BPV[t,y,s] >= -BPVMax[y,s]
          BW[t,y,s] >=  -BWMax[y,s]
          BW[t,y,s] == BWd[t,y,s] - BWc[t,y,s]
          BWe[t,y,s] <= BWMax[y,s]*DurStor

          BPV[t,y,s] == BPVd[t,y,s] - BPVc[t,y,s]
          BPVe[t,y,s] <= BPVMax[y,s]*DurStor

        end)
    end
  end
end

#Startups, shutdowns, battery's state of charge, ramp rate constraints
for y = 1:nY
  for s = 1:nS
    @constraints(m,begin
      Xc1[1,y,s] == Ec1[1,y,s] - Ac1[1,y,s]
      BWe[1,y,s] == -BWd[1,y,s] + BWc[1,y,s]*EficiencyStor
      BPVe[1,y,s] == -BPVd[1,y,s] + BPVc[1,y,s]*EficiencyStor
    end)
  end
end


for t = 2:nT
  for y = 1:nY
    for s = 1:nS
        @constraints(m,begin
          Xc1[t,y,s] == Xc1[t-1,y,s] + Ec1[t,y,s] - Ac1[t,y,s]

          Pc[t,y,s] - Pc[t-1,y,s] <= Xc1[t-1,y,s]*RampCoal*CoalCap + Ec1[t,y,s]*MinTecCoal*CoalCap
          Pc[t-1,y,s] - Pc[t,y,s] <= Xc1[t,y,s]*RampCoal*CoalCap + Ac1[t,y,s]*MinTecCoal*CoalCap

          Pg[t,y,s] - Pg[t-1,y,s] <= RampGas*GasCap
          Pg[t-1,y,s] - Pg[t,y,s] <= RampGas*GasCap

          Po[t,y,s] - Po[t-1,y,s] <= RampOil*OilCap
          Po[t-1,y,s] - Po[t,y,s] <= RampOil*OilCap

          BWe[t,y,s] == BWe[t-1,y,s] - BWd[t,y,s] + BWc[t,y,s]*EficiencyStor
          BPVe[t,y,s] == BPVe[t-1,y,s] - BPVd[t,y,s] + BPVc[t,y,s]*EficiencyStor
        end)
    end
  end
end



for t = 4:24 for y = 1:3 for s = 1:4
  @constraints(m,begin
        Xc1[t,y,s] >= sum(Ec1[k,y,s] for k = t-3:t)
        1-Xc1[t,y,s] >= sum(Ac1[k,y,s] for k = t-3:t)
      end)
end end end

for y = 1:3
  for s = 1:4
    @constraints(m,begin
      Xc1[1,y,s] >= Ec1[1,y,s]
      1-Xc1[1,y,s] >= Ac1[1,y,s]

      Xc1[2,y,s] >= sum(Ec1[k,y,s] for k = 1:2)
      1-Xc1[2,y,s] >= sum(Ac1[k,y,s] for k= 1:2)

      Xc1[3,y,s] >= sum(Ec1[k,y,s] for k = 1:3)
      1 - Xc1[3,y,s] >= sum(Ac1[k,y,s] for k=1:3)
      end)
end
end

for t = 1:nT for y = 1:1 for s = 2:nS
  @constraints(m,begin
      Pc[t,y,s] == Pc[t,y,1]
      Pc2[t,y,s] == Pc2[t,y,1]
      Pg[t,y,s] == Pg[t,y,1]
      Po[t,y,s] == Po[t,y,1]
      Xc1[t,y,s] == Xc1[t,y,1]
      Ec1[t,y,s] == Ec1[t,y,1]
      Ac1[t,y,s] == Ac1[t,y,1]
      W[t,y,s] == W[t,y,1]
      PV[t,y,s] == PV[t,y,1]
      BW[t,y,s] == BW[t,y,1]
      BPV[t,y,s] == BPV[t,y,1]
      BWc[t,y,s] == BWc[t,y,1]
      BPVc[t,y,s] == BPVc[t,y,1]
      BWd[t,y,s] == BWd[t,y,1]
      BPVd[t,y,s] == BPVd[t,y,1]
      BWe[t,y,s] == BWe[t,y,1]
      BPVe[t,y,s] == BPVe[t,y,1]
      F12[t,y,s] == F12[t,y,1]
      F23[t,y,s] == F23[t,y,1]
      F23n[t,y,s] == F23n[t,y,1]
      F13[t,y,s] == F13[t,y,1]
      Ps[t,y,s] == Ps[t,y,1]
      T1[t,y,s] == T1[t,y,1]
      T2[t,y,s] == T2[t,y,1]
      T3[t,y,s] == T3[t,y,1]
    end)
end end end

for y = 1:1 for s = 2:nS
            @constraints(m, begin
              WMax[y,s] == WMax[y,1]
              PVMax[y,s] == PVMax[y,1]
              BWMax[y,s] == BWMax[y,1]
              BPVMax[y,s] == BPVMax[y,1]
              XPS[y,s] == XPS[y,1]
            end)
end end


for t = 1:nT for y = 2:2 for s = 2:2
          @constraints(m, begin
            Pc[t,y,s] == Pc[t,y,1]
            Pc2[t,y,s] == Pc2[t,y,1]
            Pg[t,y,s] == Pg[t,y,1]
            Po[t,y,s] == Po[t,y,1]
            Xc1[t,y,s] == Xc1[t,y,1]
            Ec1[t,y,s] == Ec1[t,y,1]
            Ac1[t,y,s] == Ac1[t,y,1]
            W[t,y,s] == W[t,y,1]
            PV[t,y,s] == PV[t,y,1]
            BW[t,y,s] == BW[t,y,1]
            BPV[t,y,s] == BPV[t,y,1]
            BWc[t,y,s] == BWc[t,y,1]
            BPVc[t,y,s] == BPVc[t,y,1]
            BWd[t,y,s] == BWd[t,y,1]
            BPVd[t,y,s] == BPVd[t,y,1]
            BWe[t,y,s] == BWe[t,y,1]
            BPVe[t,y,s] == BPVe[t,y,1]
            F12[t,y,s] == F12[t,y,1]
            F23[t,y,s] == F23[t,y,1]
            F23n[t,y,s] == F23n[t,y,1]
            F13[t,y,s] == F13[t,y,1]
            Ps[t,y,s] == Ps[t,y,1]
            T1[t,y,s] == T1[t,y,1]
            T2[t,y,s] == T2[t,y,1]
            T3[t,y,s] == T3[t,y,1]
          end)
end end end

for y = 2:2 for s = 2:2
      @constraints(m,begin
      WMax[y,s] == WMax[y,1]
      PVMax[y,s] == PVMax[y,1]
      BWMax[y,s] == BWMax[y,1]
      BPVMax[y,s] == BPVMax[y,1]
      XPS[y,s] == XPS[y,1]


      WMax[y,s] >= WMax[1,1]
      PVMax[y,s] >= PVMax[1,1]
      BWMax[y,s] >= BWMax[1,1]
      BPVMax[y,s] >= BPVMax[1,1]
      XPS[y,s] >= XPS[1,1]
    end)
end end

for t = 1:nT for y = 2:2 for s = 3:3
      @constraints(m, begin
      Pc[t,y,s] == Pc[t,y,4]
      Pc2[t,y,s] == Pc2[t,y,4]
      Pg[t,y,s] == Pg[t,y,4]
      Po[t,y,s] == Po[t,y,4]
      Xc1[t,y,s] == Xc1[t,y,4]
      Ec1[t,y,s] == Ec1[t,y,4]
      Ac1[t,y,s] == Ac1[t,y,4]
      W[t,y,s] == W[t,y,4]
      PV[t,y,s] == PV[t,y,4]
      BW[t,y,s] == BW[t,y,4]
      BPV[t,y,s] == BPV[t,y,4]
      BWc[t,y,s] == BWc[t,y,4]
      BPVc[t,y,s] == BPVc[t,y,4]
      BWd[t,y,s] == BWd[t,y,4]
      BPVd[t,y,s] == BPVd[t,y,4]
      BWe[t,y,s] == BWe[t,y,4]
      BPVe[t,y,s] == BPVe[t,y,4]
      F12[t,y,s] == F12[t,y,4]
      F23[t,y,s] == F23[t,y,4]
      F23n[t,y,s] == F23n[t,y,4]
      F13[t,y,s] == F13[t,y,4]
      Ps[t,y,s] == Ps[t,y,4]
      T1[t,y,s] == T1[t,y,4]
      T2[t,y,s] == T2[t,y,4]
      T3[t,y,s] == T3[t,y,4]
      end)
end end end

for y = 2:2 for s = 3:3
      @constraints(m,begin
      WMax[y,s] == WMax[y,4]
      PVMax[y,s] == PVMax[y,4]
      BWMax[y,s] == BWMax[y,4]
      BPVMax[y,s] == BPVMax[y,4]
      XPS[y,s] == XPS[y,4]

      WMax[y,s] >= WMax[1,1]
      PVMax[y,s] >= PVMax[1,1]
      BWMax[y,s] >= BWMax[1,1]
      BPVMax[y,s] >= BPVMax[1,1]
      XPS[y,s] >= XPS[1,1]
    end)
end end

for y = 3:3 for s = 1:2
      @constraints(m,begin
      WMax[y,s] >= WMax[2,1]
      PVMax[y,s] >= PVMax[2,1]
      BWMax[y,s] >= BWMax[2,1]
      BPVMax[y,s] >= BPVMax[2,1]
      XPS[y,s] >= XPS[2,1]
    end)
end end

for y = 3:3 for s = 3:4
      @constraints(m,begin
      WMax[y,s] >= WMax[2,4]
      PVMax[y,s] >= PVMax[2,4]
      BWMax[y,s] >= BWMax[2,4]
      BPVMax[y,s] >= BPVMax[2,4]
      XPS[y,s] >= XPS[2,4]
    end)
end end

for y = 2:2 for s = 2:4
    @constraint(m, Pc2_max[y,s] == Pc2_max[y,1])
    @constraint(m, XFCo2[y,s] == XFCo2[y,1])
end end

for y = 3:3 for s = 2:2
      @constraint(m, Pc2_max[y,s] == Pc2_max[y,1])
      @constraint(m, XFCo2[y,s] == XFCo2[y,1])
end end

for y = 3:3 for s = 2:2
      @constraint(m, Pc2_max[y,s] >= Pc2_max[2,1])
      @constraint(m, XFCo2[y,s] >= XFCo2[2,1])
    end
  end

for y = 3:3 for s = 3:3
      @constraint(m, Pc2_max[y,s] == Pc2_max[y,4])
      @constraint(m, XFCo2[y,s] == XFCo2[y,4])
end end

for y = 3:3 for s = 3:3
      @constraint(m, Pc2_max[y,s] >= Pc2_max[2,1])
      @constraint(m, XFCo2[y,s] >= XFCo2[2,4])
end end


for y = 1:1 for s = 1:nS
    @constraint(m, Pc2_max[y,s] == 0)
    @constraint(m, XFCo2[y,s] == 0)
end end


for s = 1:nS
  @constraints(m, begin
    Pc2_max[2,s] >= Pc2_max[1,s]
    Pc2_max[3,s] >= Pc2_max[2,s]

    XFCo2[2,s] >= XFCo2[1,s]
    XFCo2[3,s] >= XFCo2[2,s]

    WMax[2,s] >= WMax[1,s]
    WMax[3,s] >= WMax[2,s]

    PVMax[2,s] >= PVMax[1,s]
    PVMax[3,s] >= PVMax[2,s]

    BWMax[2,s] >= BWMax[1,s]
    BWMax[3,s] >= BWMax[2,s]

    BPVMax[2,s] >= BPVMax[1,s]
    BPVMax[3,s] >= BPVMax[2,s]

    XPS[2,s] >= XPS[1,s]
    XPS[3,s] >= XPS[2,s]

#Either line or phase shifter is built
    XFCo2[3,s] + XPS[3,s] <= 1

#initial and final conditions for storage
    BWe[1,1,s] == 0
    BWe[24,1,s] == 0
    BPVe[1,1,s] == 0
    BPVe[24,1,s] == 0

    BWe[1,2,s] == 0
    BWe[24,2,s] == 0
    BPVe[1,2,s] == 0
    BPVe[24,2,s] == 0

    BWe[1,3,s] == 0
    BWe[24,3,s] == 0
    BPVe[1,3,s] == 0
    BPVe[24,3,s] == 0

  end)
end

time = @timed status = optimize!(m)

println("Objective Value = ", getobjectivevalue(m))
