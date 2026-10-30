
if nsimul0 == 0,
   vecmseinput = zeros(1,nsimul);
   vecmseitdepHTminmse = zeros(1,nsimul);
   vecmseitdepHTmingcv = zeros(1,nsimul);
   vecmsefixedHT = zeros(1,nsimul);
   vecmseitdepSTminmse = zeros(1,nsimul);
   vecmseitdepSTmingcv = zeros(1,nsimul);
   vecmsefixedST = zeros(1,nsimul);

   vecfalsepositdepHTminmse = zeros(1,nsimul);
   vecfalsepositdepHTmingcv = zeros(1,nsimul);
   vecfalseposfixedHT = zeros(1,nsimul);
   vecfalsepositdepSTminmse = zeros(1,nsimul);
   vecfalsepositdepSTmingcv = zeros(1,nsimul);
   vecfalseposfixedST = zeros(1,nsimul);

   vecfalsenegitdepHTminmse = zeros(1,nsimul);
   vecfalsenegitdepHTmingcv = zeros(1,nsimul);
   vecfalsenegfixedHT = zeros(1,nsimul);
   vecfalsenegitdepSTminmse = zeros(1,nsimul);
   vecfalsenegitdepSTmingcv = zeros(1,nsimul);
   vecfalsenegfixedST = zeros(1,nsimul);

   vectruepositdepHTminmse = zeros(1,nsimul);
   vectruepositdepHTmingcv = zeros(1,nsimul);
   vectrueposfixedHT = zeros(1,nsimul);
   vectruepositdepSTminmse = zeros(1,nsimul);
   vectruepositdepSTmingcv = zeros(1,nsimul);
   vectrueposfixedST = zeros(1,nsimul);
else
   vecmseinput = [vecmseinput zeros(nsimul-nsimul0,1)];
   vecmseitdepHTminmse = [vecmseitdepHTminmse zeros(nsimul-nsimul0,1)];
   vecmseitdepHTmingcv = [vecmseitdepHTmingcv zeros(nsimul-nsimul0,1)];
   vecmsefixedHT = [vecmsefixedHT zeros(nsimul-nsimul0,1)];
   vecmseitdepSTminmse = [vecmseitdepSTminmse zeros(nsimul-nsimul0,1)];
   vecmseitdepSTmingcv = [vecmseitdepSTmingcv zeros(nsimul-nsimul0,1)];
   vecmsefixedST = [vecmsefixedST zeros(nsimul-nsimul0,1)];

   vecfalsepositdepHTminmse = [vecfalsepositdepHTminmse zeros(nsimul-nsimul0,1)];
   vecfalsepositdepHTmingcv = [vecfalsepositdepHTmingcv zeros(nsimul-nsimul0,1)];
   vecfalseposfixedHT = [vecfalseposfixedHT zeros(nsimul-nsimul0,1)];
   vecfalsepositdepSTminmse = [vecfalsepositdepSTminmse zeros(nsimul-nsimul0,1)];
   vecfalsepositdepSTmingcv = [vecfalsepositdepSTmingcv zeros(nsimul-nsimul0,1)];
   vecfalseposfixedST = [vecfalseposfixedST zeros(nsimul-nsimul0,1)];

   vecfalsenegitdepHTminmse = [vecfalsenegitdepHTminmse zeros(nsimul-nsimul0,1)];
   vecfalsenegitdepHTmingcv = [vecfalsenegitdepHTmingcv zeros(nsimul-nsimul0,1)];
   vecfalsenegfixedHT = [vecfalsenegfixedHT zeros(nsimul-nsimul0,1)];
   vecfalsenegitdepSTminmse = [vecfalsenegitdepSTminmse zeros(nsimul-nsimul0,1)];
   vecfalsenegitdepSTmingcv = [vecfalsenegitdepSTmingcv zeros(nsimul-nsimul0,1)];
   vecfalsenegfixedST = [vecfalsenegfixedST zeros(nsimul-nsimul0,1)];

   vectruepositdepHTminmse = [vectruepositdepHTminmse zeros(nsimul-nsimul0,1)];
   vectruepositdepHTmingcv = [vectruepositdepHTmingcv zeros(nsimul-nsimul0,1)];
   vectrueposfixedHT = [vectrueposfixedHT zeros(nsimul-nsimul0,1)];
   vectruepositdepSTminmse = [vectruepositdepSTminmse zeros(nsimul-nsimul0,1)];
   vectruepositdepSTmingcv = [vectruepositdepSTmingcv zeros(nsimul-nsimul0,1)];
   vectrueposfixedST = [vectrueposfixedST  zeros(nsimul-nsimul0,1)];
end
