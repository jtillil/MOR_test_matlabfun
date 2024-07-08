%%% Version: 19 Jan 2023
%%%
%%% par  =  <MODELNAME>_parameters(model)
%%%
%%% This function creates a structure with all parameters
%%%
%%% Input :  model      model structure containing the index structure of
%%%                     the model and its initial values 
%%%                      
%%% Output : par        parameter vector
%%%
%%% References:
%%%    + Hornberg, J et al.
%%%      'Control of MAPK signaling: from complexity to what really matters'
%%%      Oncogene, Volume 24, 2005
%%%    + Schoeberl, Birgit; Eichler-Johnson, Claudia; Gilles, Ernst Dieter;
%%%      Mueller, Gertrazd
%%%      'Computational modeling of the dynamics of the MAP 
%%%      kinase casade activated by surface and internalized 
%%%      EGF receptors'
%%%      Nature Biotechnology Volume 20, 2002
%%%    + Corresponding Model from BioModels Database:
%%%      http://www.ebi.ac.uk/biomodels-main/BIOMD0000000019
%%%
%%% 
%%% Authors: Jane Knoechel and Wilhelm Huisinga
%%%

function par = Hornberg2005EGFRsignalling_parameters(model)

%%% assign model indexing
I = model.I;

%%% initializing parameter vector
par = NaN(I.npar,1);

%%% -----------------------------------------------------------------------
%%% parameter values from Hornberg et al from Supplementary table 2 
%%% in [1/s] or [1/[M*s]]

%%% CORRECTION NOTE:
%%% Some parameter values were corrected in comparison to the Hornberg et al (2005) 
%%% paper and supplement, or if not stated there, taken from 
%%% http://models.cellml.org/exposure/48c4b41256d698c4d18aace1cb159865 or 
%%% folder hornberg_binder_bruggeman_schoeberl_heinrich_westerhoff_2005-63a4a832826e 
%%% file hornberg_model_errors.pdf
%%%
%%% The correced or complemented parameters are commented 
%%% with 'CORRECTED (see note at the top)'

par(I.k1)      = 3e+7;         % second order [1/[M*s]]
par(I.kd1)     = 3.84e-3;      % first order [1/s]
par(I.k2)      = 1.66e-5;      % second order [1/[M* s]]
par(I.kd2)     = 0.1;          % first order [1/s]
par(I.k3)      = 1;            % first order [1/s]
par(I.kd3)     = 1e-2;         % first order [1/s]
par(I.k4)      = 1.73e-7;      % second order [Receptors/s]
par(I.kd4)     = 1.66e-3;      % first order [1/s]
par(I.k5)      = 0;            % 
par(I.kd5)     = 1.46e-2;      % first order [1/s]
par(I.k6)      = 5e-4;         % first order [1/s]
par(I.kd6)     = 5e-3;         % first order [1/s]
par(I.k7)      = par(I.k6);                                
par(I.kd7)     = par(I.kd6);

par(I.k8)        = 1.66e-6;      % second order [1/[M *s]]
par(I.kd8)       = 2e-1;         % first order [1/s]
par(I.k9)      = par(I.k6);                                
par(I.kd9)     = par(I.kd6);
%par(I.k10)     = 1.4e+5;       %%% PARAMETER k10 NOT EXISTENT IN THE MODEL
par(I.k10b)    = 5.43e-2;      % second order [1/[M *s]]
%par(I.kd10b)   = 0;             %%% PARAMETER k10 NOT EXISTENT IN THE MODEL 
par(I.kd10)    = 1.10e-2;      % first order [1/s]
par(I.k11)     = par(I.k2);                                
par(I.kd11)    = par(I.kd2);
par(I.k12)     = par(I.k3);                                
par(I.kd12)    = par(I.kd3);

%par(I.k13)     = 4.28e+4;     % zero order [receptors/s]
par(I.k13)     = 2.17;        % CORRECTED (see note at the top)
par(I.kd13)    = 0;            % first order [1/s]
par(I.k14)     = par(I.k8);                                
par(I.kd14)    = par(I.kd8);
par(I.k15)     = 1e+4;         % first order [1/s]
par(I.kd15)    = 0;            % first order [1/s]
par(I.k16)     = 1.66e-5;      % second order [1/[M *s]]
%par(I.kd16)    = 0;            % first order [1/s]
par(I.kd16)    = 2.75E-01;     % CORRECTED (see note at the top, taken from matlab code identical to kd63)
par(I.k17)     = 1.66e-5;      % second order [1/[M *s]]
par(I.kd17)    = 6e-2;         % first order [1/s]
par(I.k18)     = 2.5e-5;       % second order [1/[M *s]]
%par(I.kd18)    = 3.8e+4;       % first order [1/s]
par(I.kd18)    = 1.3;          % CORRECTED (see note at the top)
par(I.k19)     = 1.66e-7;      % second order [1/[M *s]]
par(I.kd19)    = 5e-1;         % first order [1/s]
par(I.k20)     = 3.5e-6;       % second order [1/[M *s]]
par(I.kd20)    = 4e-1;         % first order [1/s]
par(I.k21)     = 3.66e-7;      % second order [1/[M *s]]
par(I.kd21)    = 2.30e-2;      % first order [1/s]        % 2.2E-5 in Schoeberl2002;
par(I.k22)     = 3.5e-5;       % second order [1/[M *s]]
par(I.kd22)    = 1e-1;         % first order [1/s]
par(I.k23)     = 6;            % first order [1/s]
par(I.kd23)    = 6e-2;         % first order [1/s]
par(I.k24)     = par(I.k16);  
par(I.kd24)    = 5.5e-1;       % first order [1/s]
par(I.k25)     = 1.66e-5;      % second order [1/[M *s]]
par(I.kd25)    = 2.14e-2;      % first order [1/s]
par(I.k26)     = par(I.k18);                               
par(I.kd26)    = par(I.kd18);
par(I.k27)     = par(I.k19);                               
par(I.kd27)    = par(I.kd19);
par(I.k28)     = 1.66E-06;     % second order [1/[M *s]]
par(I.kd28)    = 5.30E-03;     % first order [1/s]
%par(I.k29)     = 1.17E-02;     % second order [1/[M *s]]
par(I.k29)     = 1.17e-6;      % CORRECTED (see note at the top)
par(I.kd29)    = 1;            % first order [1/s]
par(I.k30)     = par(I.k20);                               
par(I.kd30)    = par(I.kd20);
par(I.k31)     = par(I.k21);                               
par(I.kd31)    = par(I.kd21);
par(I.k32)     = 4.00E-07;     % second order [1/[M *s]]
par(I.kd32)    = 1.00E-01;     % first order [1/s]
par(I.k33)     = 3.50E-05;     % second order [1/[M *s]]
par(I.kd33)    = 2.00E-01;     % first order [1/s]
par(I.k34)     = 7.50E-06;     % second order [1/[M *s]]
par(I.kd34)    = 3.00E-02;     % first order [1/s]
par(I.k35)     = 7.50E-06;     % second order [1/[M *s]]
par(I.kd35)    = 1.50E-03;     % first order [1/s]
par(I.k36)     = 5.00E-03;     % first order [1/s]
par(I.kd36)    = 0;            % first order [1/s]
par(I.k37)     = 1.50E-06;     % second order [1/[M *s]]
par(I.kd37)    = 3.00E-01;     % first order [1/s]
par(I.k38)     = par(I.k16);
par(I.kd38)    = par(I.kd24);
par(I.k39)     = par(I.k37);                               
par(I.kd39)    = par(I.kd37);
par(I.k40)     = 5.00E-05;     % second order [1/[M *s]]
par(I.kd40)    = 6.40E-02;     % first order [1/s]
par(I.k41)     = 5.00E-05;     % second order [1/[M *s]]
par(I.kd41)    = 4.29E-02;     % first order [1/s]
par(I.k42)     = 1.18E-04;     % second order [1/[M *s]]
par(I.kd42)    = 2.00E-01;     % first order [1/s]
par(I.k43)     = 0;            % second order [1/[M *s]]
par(I.kd43)    = 1.00E+00;     % first order [1/s]
%par(I.k44)     = 1.95E-01;     % second order [1/[M *s]]
par(I.k44)     = 1.95e-5;      % CORRECTED (see note at the top)
par(I.kd44)    = 1.83E-02;     % first order [1/s]
par(I.k45)     = 0;            % second order [1/[M *s]]
%par(I.kd45)    = 3.81E+04;     % first order [1/s]
par(I.kd45)    = 3.5;          % CORRECTED (see note at the top)
par(I.k46)     = par(I.k44);                               
par(I.kd46)    = 3.30E-02;     % first order [1/s]  
par(I.k47)     = 0;            % second order [1/[M *s]]
%par(I.kd47)    = 3.82E+04;     % first order [1/s] 
par(I.kd47)    = 2.9;          % CORRECTED (see note at the top)
%par(I.k48)     = 2.38E-01;     % second order [1/[M *s]]
par(I.k48)     = 2.38e-5;      % CORRECTED (see note at the top)
par(I.kd48)    = 8.00E-01;     % first order [1/s]
par(I.k49)     = 0;            % second order [1/[M *s]]
%par(I.kd49)    = 5.80E-02;     % first order [1/s]
par(I.kd49)    = 5.68e-2;      % CORRECTED (see note at the top)
par(I.k50)     = 4.50E-07;     % second order [1/[M *s]]
par(I.kd50)    = 5.00E-01;     % first order [1/s]
par(I.k51)     = par(I.k49);                               
par(I.kd51)    = par(I.kd49);
%par(I.k52)     = 8.91E-01;     % second order [1/[M *s]]
par(I.k52)     = 8.91e-5;      % CORRECTED (see note at the top)
par(I.kd52)    = 3.30E-02;     % first order [1/s]
par(I.k53)     = 0;            % second order [1/[M *s]]
par(I.kd53)    = 1.60E+01;     % first order [1/s]
par(I.k54)     = par(I.k52);                               
par(I.kd54)    = 1.83E-02;     % first order [1/s]
par(I.k55)     = 0;            % second order [1/[M *s]]
%par(I.kd55)    = 3.82E+04;     % first order [1/s]
par(I.kd55)    = 5.7;          % CORRECTED (see note at the top)
par(I.k56)     = 2.35E-05;     % second order [1/[M *s]]
par(I.kd56)    = 6.00E-01;     % first order [1/s]
par(I.k57)     = 0;            % second order [1/[M *s]]
par(I.kd57)    = 2.46E-01;     % first order [1/s]
%par(I.k58)     = 8.33E-02;     % second order [1/[M *s]]
par(I.k58)     = 8.33e-6;      % CORRECTED (see note at the top)
par(I.kd58)    = 5.00E-01;     % first order [1/s]
par(I.k59)     = par(I.k57);                               
par(I.kd59)    = par(I.kd57);
par(I.k60)     = 5.50E-03;     % first order [1/s]
par(I.kd60)    = 0;            % first order [1/s]
par(I.k61)     = 6.70E-04;     % first order [1/s]
par(I.kd61)    = 0;            % first order [1/s]
par(I.k62)  = par(I.k60);                               
par(I.kd62) = par(I.kd61);
par(I.k63)  = par(I.k16);      % first order [1/s]                      
par(I.kd63)    = 2.75E-01;     % first order [1/s]
par(I.k64)  = par(I.k17);                               
par(I.kd64) = par(I.kd17);
par(I.k65)  = par(I.k18);                               
par(I.kd65) = par(I.kd18);
par(I.k66)  = par(I.k19);                               
par(I.kd66) = par(I.kd19);
par(I.k67)  = par(I.k20);                               
par(I.kd67) = par(I.kd20);
par(I.k68)  = par(I.k21);                               
par(I.kd68) = par(I.kd21);
par(I.k69)  = par(I.k22);                               
par(I.kd69) = par(I.kd22);
par(I.k70)  = par(I.k23);                               
par(I.kd70) = par(I.kd23);
par(I.k71)  = par(I.k16);                               
par(I.kd71) = par(I.kd24);
par(I.k72)  = par(I.k25);                               
par(I.kd72) = par(I.kd25);
par(I.k73)  = par(I.k18);                               
par(I.kd73) = par(I.kd18);
par(I.k74)  = par(I.k19);                               
par(I.kd74) = par(I.kd19);
par(I.k75)  = par(I.k28);                               
par(I.kd75) = par(I.kd28);
par(I.k76)  = par(I.k29);                               
par(I.kd76) = par(I.kd29);
par(I.k77)  = par(I.k20);                               
par(I.kd77) = par(I.kd20);
par(I.k78)  = par(I.k21);                               
par(I.kd78) = par(I.kd21);
par(I.k79)  = par(I.k32);                               
par(I.kd79) = par(I.kd32);
par(I.k80)  = par(I.k34);                               
par(I.kd80) = par(I.kd34);
par(I.k81)  = par(I.k37);                               
par(I.kd81) = par(I.kd37);
par(I.k82)  = par(I.k37);                               
par(I.kd82) = par(I.kd37);
par(I.k83)  = par(I.k41);                               
par(I.kd83) = par(I.kd41);
par(I.k84)  = par(I.k42);                               
par(I.kd84) = par(I.kd42);
par(I.k85)  = par(I.k43);                               
par(I.kd85) = par(I.kd43); 
par(I.k86)  = par(I.k44);                               
par(I.kd86) = par(I.kd44);
par(I.k87)  = par(I.k45);                               
par(I.kd87) = par(I.kd45); 
par(I.k88)  = par(I.k44);                               
par(I.kd88) = par(I.kd44);
par(I.k89)  = par(I.k47);                               
par(I.kd89) = par(I.kd47); 
par(I.k90)  = par(I.k48);                               
par(I.kd90) = par(I.kd48);
par(I.k91)  = par(I.k49);                               
par(I.kd91) = par(I.kd49); 
par(I.k92)  = par(I.k50);                               
par(I.kd92) = par(I.kd50);
par(I.k93)  = par(I.k49);                               
par(I.kd93) = par(I.kd49); 
par(I.k94)  = par(I.k52);                               
par(I.kd94) = 1.83E-02;     % first order [1/s]
par(I.k95)  = par(I.k53);                               
par(I.kd95) = par(I.kd53); 
par(I.k96)  = par(I.k52);                               
par(I.kd96) = 1.83E-02;     % first order [1/s]
par(I.k97)  = par(I.k55);                               
par(I.kd97)  = par(I.kd55); 
par(I.k98)  = par(I.k56);                               
par(I.kd98) = par(I.kd56);
par(I.k99)  = par(I.k57);                               
par(I.kd99)  = par(I.kd57); 
par(I.k100)  = par(I.k58);                              
par(I.kd100) = par(I.kd58);
par(I.k101)  = par(I.k57);                              
par(I.kd101)  = par(I.kd57);
par(I.k102)  = par(I.k6);                               
par(I.kd102) = par(I.kd6);
par(I.k103)  = par(I.k6);                               
par(I.kd103) = par(I.kd6);
par(I.k104)  = par(I.k6);                               
par(I.kd104) = par(I.kd6);
par(I.k105)  = par(I.k6);                               
par(I.kd105) = par(I.kd6);
par(I.k106)  = par(I.k4);                               
par(I.kd106) = par(I.kd4);
par(I.k107)  = par(I.k5);                               
par(I.kd107) = par(I.kd5);
par(I.k108)  = par(I.k6);                               
par(I.kd108) = par(I.kd6);
par(I.k109)  = par(I.k4);                               
par(I.kd109) = par(I.kd4);
par(I.k110)  = par(I.k5);                               
par(I.kd110) = par(I.kd5);
par(I.k111)  = par(I.k6);                               
par(I.kd111) = par(I.kd6);
par(I.k112)  = par(I.k4);                               
par(I.kd112) = par(I.kd4);
par(I.k113)  = par(I.k5);                               
par(I.kd113) = par(I.kd5);
par(I.k114)  = par(I.k6);                               
par(I.kd114) = par(I.kd6);
par(I.k115)  = par(I.k4);                               
par(I.kd115) = par(I.kd4);
par(I.k116)  = par(I.k5);                               
par(I.kd116) = par(I.kd5);
par(I.k117)  = par(I.k6);                               
par(I.kd117) = par(I.kd6);
par(I.k118)  = par(I.k4);                               
par(I.kd118) = par(I.kd4);
par(I.k119)  = par(I.k5);                               
par(I.kd119) = par(I.kd5);
par(I.k120)  = par(I.k6);                               
par(I.kd120) = par(I.kd6);
par(I.k121)  = par(I.k4);                               
par(I.kd121) = par(I.kd4);
par(I.k122)  = par(I.k5);                               
par(I.kd122) = par(I.kd5);
par(I.k123)  = par(I.k6);                               
par(I.kd123) = par(I.kd6);
par(I.k124)  = par(I.k4);                               
par(I.kd124) = par(I.kd4);
par(I.k125)  = par(I.k5);                               
par(I.kd125) = par(I.kd5);
par(I.k126)  = 1.66E-07;     % second order [1/[M *s]]   
par(I.kd126) = 2;            % first order [1/s]
par(I.k127)  = 1.66E-07;     % second oder [1/[M *s]]                            
par(I.kd127) = par(I.kd126);             
par(I.k128)  = par(I.k126);                             
par(I.kd128) = par(I.kd126);
par(I.k129)  = par(I.k126);                             
par(I.kd129) = par(I.kd126);
par(I.k130)  = par(I.k126);                             
par(I.kd130) = par(I.kd126);
par(I.k131)  = par(I.k126);                             
par(I.kd131) = par(I.kd126);
par(I.k132)  = par(I.k60);                              
par(I.kd132) = par(I.kd60);
par(I.k133)  = par(I.k60);                              
par(I.kd133) = par(I.kd60);
par(I.k134)  = par(I.k60);                              
par(I.kd134) = par(I.kd60);
par(I.k135)  = par(I.k60);                              
par(I.kd135) = par(I.kd60);
par(I.k136)  = par(I.k60);                              
par(I.kd136) = par(I.kd60);
par(I.k137)  = par(I.k60);                              
par(I.kd137) = par(I.kd60);
par(I.k138)  = par(I.k60);                              
par(I.kd138) = par(I.kd60);
par(I.k139)  = par(I.k60);                              
par(I.kd139) = par(I.kd60);
par(I.k140)  = par(I.k60);                              
par(I.kd140) = par(I.kd60);
par(I.k141)  = par(I.k60);                              
par(I.kd141) = par(I.kd60);
par(I.k142)  = par(I.k60);                              
par(I.kd142) = par(I.kd60);
par(I.k143)  = 0;            % second order [1/[M *s]]  
par(I.kd143) = 1.00E-04;     % first order [1/s]
par(I.k144)  = 0;   
par(I.kd144) = 1.00E-04;
par(I.k145)  = 0;   
par(I.kd145) = 1.00E-04;
par(I.k146)  = 0;   
par(I.kd146) = 1.00E-04;
par(I.k147)  = 0;   
par(I.kd147) = 1.00E-04;
par(I.k148)  = 0;   
par(I.kd148) = 1.00E-04;


end
