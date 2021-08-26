function matrixEliminated =  backSubstitution(matrixRowEchelon, ...
    indexColPivot, rankOfMatrix)
for iPivot = rankOfMatrix:-1:1
   for iRow = (iPivot-1):-1:1
       needSubtract = matrixRowEchelon(iRow, indexColPivot(iPivot))==1;
       if needSubtract
           matrixRowEchelon(iRow,:) = mod(matrixRowEchelon(iRow,:)+ ...
               matrixRowEchelon(iPivot,:),2);
       end
   end
end
matrixEliminated = matrixRowEchelon;
end