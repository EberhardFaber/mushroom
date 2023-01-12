function varargout = untitled(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 12-Jan-2023 03:23:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% Аэродинамические параметры 
% Константы
global Fsopr VklSoprVozduha y E
global T betta epsilong psi D Rvkl_RGS theta epsilonv g m FY
global bettatr lambda V Tc deltax ABg ABv
global Mvozd Cx0 Cd R StandartAtmosph Smid
global t k mrs mt Ttyagi
Mvozd=0.0289644; % Молярная масса воздуха, [кг/моль]
g=9.81; % Ускорение свободного падения, [м/с2]
T0=288;% Температура воздуха над морем, [К]
L=-0.0065; % Температурный градиент, [К/м]
rovozd(1)=1.225; % Плотность воздуха над морем, [кг/м3]

R=8.31447; % Универсальная газовая постоянная, [Дж/(моль*К)]

StandartAtmosph=[0	500	1000	1500	2000	2500	3000	4000	5000	6000	7000	8000	9000	10000	11000	12000	14000	16000	18000	20000	24000	28000	32000	36000	40000;
288.2	284.9	281.7	278.4	275.2	271.9	268.7	262.2	255.7	249.2	242.7	236.2	229.7	223.3	216.8	216.7	216.7	216.7	216.7	216.7	220.6	224.5	228.5	239.3	250.4;
101330	95464	89877	84559	79499	74690	70123	61661	54052	47217	41106	35653	30801	26500	22700	19399	14170	10353	7565	5529	2971	1616	889	499	287;
1.225	1.1673	1.1117	1.0581	1.0065	0.9569	0.9083	0.8194	0.7365	0.6601	0.59	0.5258	0.4671	0.4135	0.3648	0.3119	0.2279	0.1665	0.1216	0.0889	0.0469	0.0251	0.0136	0.00726	0.004
];

% Параметры ракеты РВВ-АЕ
Smid = pi()/4*0.2^2;


mrs=175; % Старотвая масса ракеты (кг)
mt=58.5; % Масса топлива
Ttyagi=3.5; % Время работы двигателя
global Lhordi Lkrila Skrila Akrila

Lhordi=0.7;
Lkrila=0.53;
Skrila=Lhordi*Lkrila;
Akrila=Lkrila/Lhordi;
% Условия моделирования

% Параметры модели
% T=0.1;
KM=1000000;

% Параметры РГС

% k_vkl_AN=ceil(Ttyagi/T);
deltax=5;
lambda=0.008;
Tc=0.2;


% [x,y,z, xc, yc ,zc, D]=deal(zeros(1,KM));

[t, mrt]=deal(zeros(1,KM));

%Параметры цели
[xc, yc ,zc, Vcz]=deal(zeros(1,KM));
xc(1)=str2double(get(handles.edit3, 'String'));
yc(1)=str2double(get(handles.edit4, 'String'));
zc(1)=str2double(get(handles.edit5, 'String'));
psic=0*pi()/180;
thetac=0*pi()/180;
Vc=0;

[x, y ,z, E, test, D, psi, theta, V]=deal(zeros(1,KM));
% Условия пуска
% Координаты (м)
initial_rocket_y=str2double(get(handles.edit2, 'String'));
x(1)=0;
y(1)=initial_rocket_y;
z(1)=0;
V(1)=70; % Модуль начальной скорости (м/с)
E(1)=mrs*g*y(1)+(mrs*V(1)^2)/2; % Полная энергия

psi(1)=-0*pi()/180;
theta(1)=0*pi()/180;

Rakurs=1;
Rvkl_RGS=10000;
Rvikl_RGS=500;
Rvikl_UPR=50;

VklSoprVozduha=1;
VklUprV=1;
VklUprG=1;

Cx0=0.3;
Cd = [0.3 0.3 0.3 0.5 1.1 1.3 1.3 1.15 0.95 0.75 0.6 0.5 0.45 0.4 0.37 0.35;
0 0.4 0.8 0.9 0.95 1 1.1 1.15 1.2 1.4 1.6 1.8 2	2.5	3 4];

[epsilonv, epsilong, ABg, ABv, betta, bettatr, m, FY, Fsopr]=deal(zeros(1,KM));
m(1)=mrs;


k_vkl=-1;
k=1;

abg_limit=20; % Предел перегрузки по горизонтали по модулю
abv_limit=20; % Предел перегрузки по вертикали по модулю

yc(1)=0;

msize=10; % Размерность матрицы промахов

xnach=10000;
xkon=35000;
xstep=(xkon-xnach)/msize;

znach=-25000;
zkon=25000;
zstep=(zkon-znach)/msize;

initial_T=0.01;
secondary_T=0.0001; 
T=initial_T;

promah=zeros((xkon-xnach)/xstep, (zkon-znach)/zstep);
square_area=((xkon-xnach)/msize)*((zkon-znach)/msize); % Площадь одного "квадрата"
distruction_area=0; % Общая площадь поражения
allowable_miss=1; % Максимально допустимый промах
step_change_distance=50; %Расстояние до цели, на котором меняется шаг интегрирования



% for i=xnach:xstep:xkon-0.01 % Координата цели по x
%     xc(1)=i+xstep/2;
%     for j=znach:zstep:zkon-0.01 % Координата цели по z
%         zc(1)=j+zstep/2;
        while y(k)>0
            k=k+1;
            t(k)=(k-1)*T;
            m(k)=mp(t(k-1));
            D(k-1)=sqrt((xc(k-1)-x(k-1))^2+(yc(k-1)-y(k-1))^2+(zc(k-1)-z(k-1))^2);
            if D(k-1)>step_change_distance
                T=initial_T;
            else
                T=secondary_T;
            end
            if k>2
                if D(k-1)>D(k-2)
                    break
                end
            end
            if k_vkl<0 && D(k-1)<=Rvkl_RGS
               k_vkl=k-1;
            end
            epsilonv(k-1)=asin((yc(k-1)-y(k-1))/D(k-1));
            epsilong(k-1)=atan2(zc(k-1)-z(k-1),xc(k-1)-x(k-1));
            if (k-1)*T<Ttyagi %Если двигатель работает
                abg=0;
            else
                betta(k-1)=epsilong(k-1)-psi(k-1);
                if D(k-1)>Rvkl_RGS % Наведение на автономном участке
                    bettatr(k-1)=asin((lambda*D(k-1))/(2*V(k-1)*Tc*deltax));
                    if zc(1)-z(1)>=0
                        abg=V(k-1)/0.5*(betta(k-1)-bettatr(k-1));
                    else
                        abg=V(k-1)/0.5*(betta(k-1)+bettatr(k-1));
                    end
                else %Самонаведение
                    bettatr(k-1)=asin((lambda*D(k-1))/(2*V(k-1)*Tc*deltax));
                    if zc(1)-z(1)>=0
                        abg=V(k-1)/(2*T)*(betta(k-1)-bettatr(k-1))*2;
                    else
                        abg=V(k-1)/(2*T)*(betta(k-1)+bettatr(k-1))*2;
                    end
                end
            end
            if abs(abg)>abg_limit
                abg=abg/(abs(abg)/abg_limit);
            end
            ABg(k-1)=VklUprG*abg; % Управление по горизонтали
            if (k-1)*T<Ttyagi
                abv=g*cos(theta(k-1));
            else
                abv=(epsilonv(k-1)-epsilonv(k-2))/T*10000+g*cos(theta(k-1));
            end
            if abs(abg)>abv_limit
                abg=abg/(abs(abg)/abv_limit);
            end
            ABv(k-1)=VklUprV*abv; %Управление по вертикали

            xc(k)=xc(k-1)+Vc*T*cos(psic)*cos(thetac);
            yc(k)=yc(k-1)+Vc*T*sin(thetac);
            zc(k)=zc(k-1)+Vc*T*sin(psic)*cos(thetac);

            x(k)=x(k-1)+V(k-1)*T*cos(psi(k-1))*cos(theta(k-1));
            y(k)=y(k-1)+V(k-1)*T*sin(theta(k-1));
            z(k)=z(k-1)+V(k-1)*T*sin(psi(k-1))*cos(theta(k-1));  
            psi(k)=psi(k-1)+ABg(k-1)/(V(k-1)*cos(theta(k-1)))*T;
            theta(k)=theta(k-1)+(ABv(k-1)-g*cos(theta(k-1)))/V(k-1)*T;
            
            if imag(y(k-1))~=0 || imag(V(k-1))~=0
               fprintf('warning! комплексные числа! \n')
               break
            end
            V(k)=v();
                    
        end
        

%         if D(k-1)>allowable_miss
%             promah((xkon-xnach)/xstep-abs(xnach-i)/xstep,1+abs(znach-j)/zstep)=allowable_miss;
%         else
%             promah((xkon-xnach)/xstep-abs(xnach-i)/xstep,1+abs(znach-j)/zstep)=D(k-1);
%             distruction_area=distruction_area+square_area;
%         end
%         promah((xkon-xnach)/xstep-abs(xnach-i)/xstep,msize-abs(znach-j)/zstep)=promah((xkon-xnach)/xstep-abs(xnach-i)/xstep,1+abs(znach-j)/zstep);
%         promah
%         fprintf('zc=%f; xc=%f; Промах (м): %.5f \n Общая площадь поражения (м^2): %.5f \n', zc(1),xc(1),D(k-1), distruction_area);
%         plot3(z(1:k-3),x(1:k-3),y(1:k-3))
%         hold on
%         plot3(-z(1:k-3),x(1:k-3),y(1:k-3))
%         plot3(zc(1),xc(1),yc(1),'r*')
%         plot3(-zc(1),xc(1),yc(1),'r*')
        
        plot(z(1:k-3),x(1:k-3))
        hold on
%         plot(-z(1:k-3),x(1:k-3))
        plot(zc(1),xc(1),'r*')
%         plot(-zc(1),xc(1),'r*')
        
%         axis equal
%         xlim([znach-10 zkon+10])
%         ylim([xnach-15 xkon+5])

%     end
% end

t = t(1:k);
V = V(1:k);
x=x(1:k);
y=y(1:k);
z=z(1:k);
xc=xc(1:k);
yc=yc(1:k);
zc=zc(1:k);

ABg=ABg(1:k);
ABv=ABv(1:k);
betta=betta(1:k);
bettatr=bettatr(1:k);
D=D(1:k);
E=E(1:k);
epsilong=epsilong(1:k);
epsilonv=epsilonv(1:k);
Fsopr=Fsopr(1:k);
FY=FY(1:k);
m=m(1:k);
psi=psi(1:k);
theta=theta(1:k);

string1 = sprintf('%.5f',D(k-1));
set(handles.edit6, 'String', string1);





function Mp = mp(tt)
global Ttyagi mrs mt
if tt<=Ttyagi
    Mp=mrs-mt*tt/Ttyagi;
else
    Mp=mrs-mt;
end


function ftyagi = Ftyagi(t)
global Ttyagi
if t<=Ttyagi
    ftyagi=40000;
else
    ftyagi=0;
end


function mach = Mach(h,v)
global R Mvozd
mach=v./sqrt(1.4*R*Temp(h)./Mvozd);


function temp = Temp(h)
global StandartAtmosph
temp=interp1(StandartAtmosph(1,1:length(StandartAtmosph)), StandartAtmosph(2,1:length(StandartAtmosph)), h);


function cx = Cx(M)
global Cd
cx=interp1(Cd(2,1:length(Cd)), Cd(1,1:length(Cd)), M);


function Ro_vozd = ro_vozd(h)
global StandartAtmosph
Ro_vozd=interp1(StandartAtmosph(1,1:length(StandartAtmosph)), StandartAtmosph(4,1:length(StandartAtmosph)), h);


function FX1l = FX1L(v, h)
global Smid
FX1l = Cx(Mach(h,v)).*(ro_vozd(h)/2.*v.^2).*Smid;


% function FX0l = FX0L(v, h)
% global Smid Cx0
% FX0l = Cx0.*(ro_vozd(h)/2.*v.^2).*Smid;
% end

function vv=v()
global m k T FY ABv ABg Fsopr VklSoprVozduha g V y E
m(k)=mp((k)*T);
FY(k-1)=m(k-1)*sqrt((ABv(k-1))^2+(ABg(k-1))^2);
Fsopr(k-1)=VklSoprVozduha*(FX1L(V(k-1),y(k-1))+FXi(V(k-1), y(k-1), FY(k-1)));
E(k)=m(k-1)*g*y(k-1)+e();
vv=sqrt(2/m(k)*(E(k)-m(k)*g*y(k)));


function ee=e()
global V k Fsopr T m
vvv=V(k-1)+(Ftyagi((k)*T)-Fsopr(k-1))/m(k-1)*T;
ee=(m(k-1)*vvv^2)/2;


function fxi = FXi(v, h, Y)
global Skrila Akrila
fxi=1/(pi()*Akrila)*Y^2/(ro_vozd(h)*v^2/2*Skrila);




% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
cla reset
string1 = sprintf('');
set(handles.edit6, 'String', string1);
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
