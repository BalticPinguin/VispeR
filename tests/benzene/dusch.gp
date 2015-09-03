set yrange [-0.1:4.4]
plot "DR.spe" u 1:2 w l title "OPA", "DR_full.spe" u 1:2 w l title "full", "changed_full.spe" u 1:2 w l title "FC_full", "changed.spe" u 1:2 w l title "FC_OPA"
pause -1
