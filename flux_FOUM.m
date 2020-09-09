%% C�lculo do fluxo na face por upwind de primeira ordem
function [dfw_vols] = flux_FOUM(M, fw)
    %% Obten��o do fluxo na face de acordo com a dire��o do fluxo
    fw_face_neig = fw(M.faces_conec);

    fw_face = zeros(length(fw_face_neig(:,1)),1);
    fw_face(fw_face_neig(:,2)>=fw_face_neig(:,1)) = fw_face_neig(fw_face_neig(:,2)>=fw_face_neig(:,1),2);
    fw_face(fw_face_neig(:,2)<fw_face_neig(:,1)) = fw_face_neig(fw_face_neig(:,2)<fw_face_neig(:,1),1);
    
    % Abaixo � mostrado um artificio usando matriz esparsa pra calcular o
    % fluxo fracional no volume (dfw_vols) como o balan�o do que entra 
    % menos o que sai pelas faces. 
    % Esse artificio foi utilizado aqui considerando que podem 
    % haver problemas 2-D e 3-D tamb�m.
    % Esse artificio tamb�m � usado pra calcular o delta de
    % satura��o (dSw_vols) no volume.
    conec = M.faces_conec;

    cph = 1;

    lines = cat(1,repmat(cph,length(conec(:,1)),1), repmat(cph, length(conec(:,2)),1));
    cols = cat(1,conec(:,1),conec(:,2));
    data = cat(1,-fw_face, fw_face);
    dfw_vols = sparse(lines,cols,data);
end