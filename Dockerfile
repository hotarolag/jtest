FROM julia:1.0.0

RUN julia -e "using Pkg; Pkg.add(\"JuMP\"); using JuMP"
RUN julia -e "using Pkg; Pkg.add(\"CSV\"); using CSV"
RUN julia -e "using Pkg; Pkg.add(\"Cbc\"); using Cbc"

RUN mkdir -p /home/script/

COPY juliatest /home/script/
ENTRYPOINT ["julia", "/home/script/SmartgridPlanning.jl"]
