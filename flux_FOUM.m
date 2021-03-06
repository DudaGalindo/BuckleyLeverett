function [dfwdSw, dSw_vols] = flux_FOUM(M, fw, Sw)
    fw_face_neig = fw(M.faces_conec);
    Sw_face_neig = Sw(M.faces_conec);
    
    fw_face = zeros(length(fw_face_neig(:,1)),1);
    fw_face(fw_face_neig(:,2)>=fw_face_neig(:,1)) = fw_face_neig(fw_face_neig(:,2)>=fw_face_neig(:,1),2);
    fw_face(fw_face_neig(:,2)<fw_face_neig(:,1)) = fw_face_neig(fw_face_neig(:,2)<fw_face_neig(:,1),1);
    
    Sw_face(fw_face_neig(:,2)>=fw_face_neig(:,1)) = Sw_face_neig(fw_face_neig(:,2)>=fw_face_neig(:,1),2);
    Sw_face(fw_face_neig(:,2)<fw_face_neig(:,1)) = Sw_face_neig(fw_face_neig(:,2)<fw_face_neig(:,1),1);
    
    conec = M.faces_conec;
    
    % Abaixo � mostrado um artificio usando matriz esparsa pra calcular o
    % fluxo fracional no volume como balan�o do que entra menos o que sai
    % pelas faces. Esse artificio tamb�m � usado pra calcular o delta de
    % satura��o (dSw) no volume.
    cph = 1;
    lines = cat(1,repmat(cph,length(conec(:,1)),1), repmat(cph, length(conec(:,2)),1));
    cols = cat(1,conec(:,1),conec(:,2));
    data = cat(1,-fw_face, fw_face);
    dfw_vols = sparse(lines,cols,data);
    
    lines = cat(1,repmat(cph,length(conec(:,1)),1), repmat(cph, length(conec(:,2)),1));
    cols = cat(1,conec(:,1),conec(:,2));
    data = cat(1,-Sw_face, Sw_face);
    dSw_vols = sparse(lines,cols,data);
    
    dfwdSw = dfw_vols./dSw_vols;
    dfwdSw(dSw_vols==0)=0;
end