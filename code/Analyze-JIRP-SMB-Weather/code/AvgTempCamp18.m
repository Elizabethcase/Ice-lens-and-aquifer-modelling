%finds average temeperatures at Camp 18 (closest to the divide)
%uses lapse rate to convert temps to divide elevation according to McNiel
%et al 2019

%%lapse rate conversion
% ADD

%% daily
C18daily = readtable('/Users/elizabeth/Documents/home_research/projects/JIRP Firn Aquifer/data/juneauIceField_weather_v1.0/LVL2/juneauicefieldCamp18AWS_daily_LVL2.csv');

figure(1)
plot(C18daily.Date,C18daily.site_temp_min,'b');
hold on;
plot(C18daily.Date,C18daily.site_temp_max,'r');

avgDailyTC18 = (mean(C18daily.site_temp_max,'omitnan') + mean(C18daily.site_temp_min,'omitnan'))/2;

%% hourly
C18hourly = readtable('/Users/elizabeth/Documents/home_research/projects/JIRP Firn Aquifer/data/juneauIceField_weather_v1.0/LVL2/juneauicefieldCamp18AWS_hourly_LVL2.csv');
C18hourly.Date = datetime(C18hourly.local_time,'InputFormat','yyyy/MM/dd HH:mm')
figure(2)
plot(C18hourly.Date,C18hourly.site_temp_min,'b');
hold on;
plot(C18hourly.Date,C18hourly.site_temp_max,'r');

avgHourlyTC18 = (mean(C18hourly.site_temp_max,'omitnan') + mean(C18hourly.site_temp_min,'omitnan'))/2;