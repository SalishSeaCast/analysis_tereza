clear

load tjarnikova-Jarnikova_Canadian_Arctic_DMS_supldata-2455bdc/osscar_mims_data.mat

coresp_mims = zeros(1,344);

mimsdate = mims.mdate
mimsdms = mims.dms
for i = 1:344
    osdate = osscar_DMS.mdate(i);
    w = mimsdms(abs(mimsdate-osdate) <0.0055);
    coresp_mims(i) = mean(w)
end

%plot(osscar_DMS.UW_corrected,mimsdms)

indic = []
tval = []
for i = 1:344
    osdms = osscar_DMS.UW_corrected(i);
    mimsdms = coresp_mims(i);
    if abs((osdms-mimsdms)/osdms) < 0.4
        tval(end+1) = (abs((osdms-mimsdms)/osdms))
        indic(end+1) = i
    end
end

osdms = osscar_DMS.UW_corrected
mimsdms = coresp_mims

w = [1:1:344]
figure
hold on
plot(w,osdms)
plot(w,mimsdms)