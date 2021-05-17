function [JB,GBE,GBC,JC,GCE,GCC,CBE,CBC]=BJT(UBE,UBC,IS,NF,NR,UT,BF,BR,TF,TR,CBE0,CBC0,U0,m,typ)
%paramjetry funkcji: typ = 1 (NPPN) -1 (PNP)
%                    UBE, UBC  napiecia baza-emiter, baza-kolektor
%                    IS prad nasycenia
%                    NF nieidealnoœæ nachylenia w przod
%                    NR nieidealnoœæ nachylenia wstecz
%                    UT potencjal elektrotermiczny
%                    BF wzmocnienie pradowe w przód
%                    BR wzmocnienie pr¹dowe wstecz
%                    TT czas przelotu w przód
%                    TR czas przelotu wstecz
%                    CBE0, CBC0 pojemnoœæ BE i BC przy zerowych napieciach
%                    U0
%                    m
UBE=UBE*typ; UBC=UBC*typ;
IB=(IS/BF)*(expo(UBE/NF/UT)-1)+(IS/BR)*(expo(UBC/NR/UT)-1);
IC=IS*(expo(UBE/NF/UT)-expo(UBC/NR/UT))-(IS/BR)*(expo(UBC/NR/UT));
GBE=(IS/BF/NF/UT)*dxpo(UBE/NF/UT);
GBC=(IS/BR/NR/UT)*dxpo(UBC/NR/UT);
GCE=(IS/NF/UT)*dxpo(UBE/NF/UT);
GCC=-(IS/NR/UT)*(1+1/BR)*dxpo(UBC/NR/UT);
JB=typ*(IB-GBE*UBE-GBC*UBC);
JC=typ*(IC-GCE*UBE-GCC*UBC);
CBE=(TF*IS/(NF*UT))*dxpo(UBE/(NF*UT))+cpx1am(UBE,CBE0,0.6,m);
CBC=(TR*IS/(NR*UT))*dxpo(UBC/(NR*UT))+cpx1am(UBC,CBC0,U0,m);