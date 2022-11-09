import itertools

import pytest
import timeout_decorator

from mavis.annotate.base import ReferenceName
from mavis.annotate.protein import Domain, DomainRegion, calculate_orf
from mavis.annotate.variant import IndelCall

from ...util import get_data


class TestDomainAlignSeq:
    def test_large_combinations_finishes_with_error(self):
        input_seq = (
            'MADDEDYEEVVEYYTEEVVYEEVPGETITKIYETTTTRTSDYEQSETSKPALAQPALAQPASAKPVERRKVIRKKVDPSK'
            'FMTPYIAHSQKMQDLFSPNKYKEKFEKTKGQPYASTTDTPELRRIKKVQDQLSEVKYRMDGDVAKTICHVDEKAKDIEHA'
            'KKVSQQVSKVLYKQNWEDTKDKYLLPPDAPELVQAVKNTAMFSKKLYTEDWEADKSLFYPYNDSPELRRVAQAQKALSDV'
            'AYKKGLAEQQAQFTPLADPPDIEFAKKVTNQVSKQKYKEDYENKIKGKWSETPCFEVANARMNADNISTRKYQEDFENMK'
            'DQIYFMQTETPEYKMNKKAGVAASKVKYKEDYEKNKGKADYNVLPASENPQLRQLKAAGDALSDKLYKENYEKTKAKSIN'
            'YCETPKFKLDTVLQNFSSDKKYKDSYLKDILGHYVGSFEDPYHSHCMKVTAQNSDKNYKAEYEEDRGKGFFPQTITQEYE'
            'AIKKLDQCKDHTYKVHPDKTKFTQVTDSPVLLQAQVNSKQLSDLNYKAKHESEKFKCHIPPDTPAFIQHKVNAYNLSDNL'
            'YKQDWEKSKAKKFDIKVDAIPLLAAKANTKNTSDVMYKKDYEKNKGKMIGVLSINDDPKMLHSLKVAKNQSDRLYKENYE'
            'KTKAKSMNYCETPKYQLDTQLKNFSEARYKDLYVKDVLGHYVGSMEDPYHTHCMKVAAQNSDKSYKAEYEEDKGKCYFPQ'
            'TITQEYEAIKKLDQCKDHTYKVHPDKTKFTAVTDSPVLLQAQLNTKQLSDLNYKAKHEGEKFKCHIPADAPQFIQHRVNA'
            'YNLSDNVYKQDWEKSKAKKFDIKVDAIPLLAAKANTKNTSDVMYKKDYEKSKGKMIGALSINDDPKMLHSLKTAKNQSDR'
            'EYRKDYEKSKTIYTAPLDMLQVTQAKKSQAIASDVDYKHILHSYSYPPDSINVDLAKKAYALQSDVEYKADYNSWMKGCG'
            'WVPFGSLEMEKAKRASDILNEKKYRQHPDTLKFTSIEDAPITVQSKINQAQRSDIAYKAKGEEIIHKYNLPPDLPQFIQA'
            'KVNAYNISENMYKADLKDLSKKGYDLRTDAIPIRAAKAARQAASDVQYKKDYEKAKGKMVGFQSLQDDPKLVHYMNVAKI'
            'QSDREYKKDYEKTKSKYNTPHDMFNVVAAKKAQDVVSNVNYKHSLHHYTYLPDAMDLELSKNMMQIQSDNVYKEDYNNWM'
            'KGIGWIPIGSLDVEKVKKAGDALNEKKYRQHPDTLKFTSIVDSPVMVQAKQNTKQVSDILYKAKGEDVKHKYTMSPDLPQ'
            'FLQAKCNAYNISDVCYKRDWYDLIAKGNNVLGDAIPITAAKASRNIASDYKYKEAYEKSKGKHVGFRSLQDDPKLVHYMN'
            'VAKLQSDREYKKNYENTKTSYHTPGDMVSITAAKMAQDVATNVNYKQPLHHYTYLPDAMSLEHTRNVNQIQSDNVYKDEY'
            'NSFLKGIGWIPIGSLEVEKVKKAGDALNERKYRQHPDTVKFTSVPDSMGMVLAQHNTKQLSDLNYKVEGEKLKHKYTIDP'
            'ELPQFIQAKVNALNMSDAHYKADWKKTIAKGYDLRPDAIPIVAAKSSRNIASDCKYKEAYEKAKGKQVGFLSLQDDPKLV'
            'HYMNVAKIQSDREYKKGYEASKTKYHTPLDMVSVTAAKKSQEVATNANYRQSYHHYTLLPDALNVEHSRNAMQIQSDNLY'
            'KSDFTNWMKGIGWVPIESLEVEKAKKAGEILSEKKYRQHPEKLKFTYAMDTMEQALNKSNKLNMDKRLYTEKWNKDKTTI'
            'HVMPDTPDILLSRVNQITMSDKLYKAGWEEEKKKGYDLRPDAIAIKAARASRDIASDYKYKKAYEQAKGKHIGFRSLEDD'
            'PKLVHFMQVAKMQSDREYKKGYEKSKTSFHTPVDMLSVVAAKKSQEVATNANYRNVIHTYNMLPDAMSFELAKNMMQIQS'
            'DNQYKADYADFMKGIGWLPLGSLEAEKNKKAMEIISEKKYRQHPDTLKYSTLMDSMNMVLAQNNAKIMNEHLYKQAWEAD'
            'KTKVHIMPDIPQIILAKANAINMSDKLYKLSLEESKKKGYDLRPDAIPIKAAKASRDIASDYKYKYNYEKGKGKMVGFRS'
            'LEDDPKLVHSMQVAKMQSDREYKKNYENTKTSYHTPADMLSVTAAKDAQANITNTNYKHLIHKYILLPDAMNIELTRNMN'
            'RIQSDNEYKQDYNEWYKGLGWSPAGSLEVEKAKKATEYASDQKYRQHPSNFQFKKLTDSMDMVLAKQNAHTMNKHLYTID'
            'WNKDKTKIHVMPDTPDILQAKQNQTLYSQKLYKLGWEEALKKGYDLPVDAISVQLAKASRDIASDYKYKQGYRKQLGHHV'
            'GFRSLQDDPKLVLSMNVAKMQSEREYKKDFEKWKTKFSSPVDMLGVVLAKKCQELVSDVDYKNYLHQWTCLPDQNDVVQA'
            'KKVYELQSENLYKSDLEWLRGIGWSPLGSLEAEKNKRASEIISEKKYRQPPDRNKFTSIPDAMDIVLAKTNAKNRSDRLY'
            'REAWDKDKTQIHIMPDTPDIVLAKANLINTSDKLYRMGYEELKRKGYDLPVDAIPIKAAKASREIASEYKYKEGFRKQLG'
            'HHIGARNIEDDPKMMWSMHVAKIQSDREYKKDFEKWKTKFSSPVDMLGVVLAKKCQTLVSDVDYKNYLHQWTCLPDQSDV'
            'IHARQAYDLQSDNLYKSDLQWLKGIGWMTSGSLEDEKNKRATQILSDHVYRQHPDQFKFSSLMDSIPMVLAKNNAITMNH'
            'RLYTEAWDKDKTTVHIMPDTPEVLLAKQNKVNYSEKLYKLGLEEAKRKGYDMRVDAIPIKAAKASRDIASEFKYKEGYRK'
            'QLGHHIGARAIRDDPKMMWSMHVAKIQSDREYKKDFEKWKTKFSSPVDMLGVVLAKKCQTLVSDVDYKNYLHQWTCLPDQ'
            'SDVIHARQAYDLQSDNMYKSDLQWMRGIGWVSIGSLDVEKCKRATEILSDKIYRQPPDRFKFTSVTDSLEQVLAKNNAIT'
            'MNKRLYTEAWDKDKTQIHIMPDTPEIMLARMNKINYSESLYKLANEEAKKKGYDLRSDAIPIVAAKASRDIISDYKYKDG'
            'YCKQLGHHIGARNIEDDPKMMWSMHVAKIQSDREYKKDFEKWKTKFSSPVDMLGVVLAKKCQTLVSDVDYKNYLHEWTCL'
            'PDQSDVIHARQAYDLQSDNIYKSDLQWLRGIGWVPIGSMDVVKCKRATEILSDNIYRQPPDKLKFTSVTDSLEQVLAKNN'
            'ALNMNKRLYTEAWDKDKTQIHIMPDTPEIMLARQNKINYSETLYKLANEEAKKKGYDLRSDAIPIVAAKASRDVISDYKY'
            'KDGYRKQLGHHIGARNIEDDPKMMWSMHVAKIQSDREYKKDFEKWKTKFSSPVDMLGVVLAKKCQTLVSDVDYKNYLHEW'
            'TCLPDQNDVIHARQAYDLQSDNIYKSDLQWLRGIGWVPIGSMDVVKCKRAAEILSDNIYRQPPDKLKFTSVTDSLEQVLA'
            'KNNALNMNKRLYTEAWDKDKTQVHIMPDTPEIMLARQNKINYSESLYRQAMEEAKKEGYDLRSDAIPIVAAKASRDIASD'
            'YKYKEAYRKQLGHHIGARAVHDDPKIMWSLHIAKVQSDREYKKDFEKYKTRYSSPVDMLGIVLAKKCQTLVSDVDYKHPL'
            'HEWICLPDQNDIIHARKAYDLQSDNLYKSDLEWMKGIGWVPIDSLEVVRAKRAGELLSDTIYRQRPETLKFTSITDTPEQ'
            'VLAKNNALNMNKRLYTEAWDNDKKTIHVMPDTPEIMLAKLNRINYSDKLYKLALEESKKEGYDLRLDAIPIQAAKASRDI'
            'ASDYKYKEGYRKQLGHHIGARNIKDDPKMMWSIHVAKIQSDREYKKEFEKWKTKFSSPVDMLGVVLAKKCQILVSDIDYK'
            'HPLHEWTCLPDQNDVIQARKAYDLQSDAIYKSDLEWLRGIGWVPIGSVEVEKVKRAGEILSDRKYRQPADQLKFTCITDT'
            'PEIVLAKNNALTMSKHLYTEAWDADKTSIHVMPDTPDILLAKSNSANISQKLYTKGWDESKMKDYDLRADAISIKSAKAS'
            'RDIASDYKYKEAYEKQKGHHIGAQSIEDDPKIMCAIHAGKIQSEREYKKEFQKWKTKFSSPVDMLSILLAKKCQTLVTDI'
            'DYRNYLHEWTCMPDQNDIIQAKKAYDLQSDSVYKADLEWLRGIGWMPEGSVEMNRVKVAQDLVNERLYRTRPEALSFTSI'
            'VDTPEVVLAKANSLQISEKLYQEAWNKDKSNITIPSDTPEMLQAHINALQISNKLYQKDWNDAKQKGYDIRADAIEIKHA'
            'KASREIASEYKYKEGYRKQLGHHMGFRTLQDDPKSVWAIHAAKIQSDREYKKAYEKSKGIHNTPLDMMSIVQAKKCQVLV'
            'SDIDYRNYLHQWTCLPDQNDVIQAKKAYDLQSDNLYKSDLEWLKGIGWLPEGSVEVMRVKNAQNLLNERLYRIKPEALKF'
            'TSIVDTPEVIQAKINAVQISEPLYRDAWEKEKANVNVPADTPLMLQSKINALQISNKRYQQAWEDVKMTGYDLRADAIGI'
            'QHAKASRDIASDYLYKTAYEKQKGHYIGCRSAKEDPKLVWAANVLKMQNDRLYKKAYNDHKAKISIPVDMVSISAAKEGQ'
            'ALASDVDYRHYLHHWSCFPDQNDVIQARKAYDLQSDSVYKADLEWLRGIGWMPEGSVEMNRVKVAQDLVNERLYRTRPEA'
            'LSFTSIVDTPEVVLAKANSLQISEKLYQEAWNKDKSNITIPSDTPEMLQAHINALQISNKLYQKDWNDTKQKGYDIRADA'
            'IEIKHAKASREIASEYKYKEGYRKQLGHHMGFRTLQDDPKSVWAIHAAKIQSDREYKKAYEKSKGIHNTPLDMMSIVQAK'
            'KCQVLVSDIDYRNYLHQWTCLPDQNDVIQAKKAYDLQSDNLYKSDLEWLKGIGWLPEGSVEVMRVKNAQNLLNERLYRIK'
            'PEALKFTSIVDTPEVIQAKINAVQISEPLYRNAWEKEKANVNVPADTPLMLQSKINALQISNKRYQQAWEDVKMTGYDLR'
            'ADAIGIQHAKASRDIASDYLYKTAYEKQKGHYIGCRSAKEDPKLVWAANVLKMQNDRLYKKAYNDHKAKISIPVDMVSIS'
            'AAKEGQALASDVDYRHYLHHWSCFPDQNDVIQARKAYDLQSDSVYKADLEWLRGIGWMPEGSVEMNRVKVAQDLVNERLY'
            'RTRPEALSFTSIVDTPEVVLAKANSLQISEKLYQEAWNKDKSNITIPSDTPEMLQAHINALQISNKLYQKDWNDTKQKGY'
            'DIRADAIEIKHAKASREIASEYKYKEGYRKQLGHHMGFRTLQDDPKSVWAIHAAKIQSDREYKKAYEKSKGIHNTPLDMM'
            'SIVQAKKCQVLVSDIDYRNYLHQWTCLPDQNDVIQAKKAYDLQSDNLYKSDLEWLKGIGWLPEGSVEVMRVKNAQNLLNE'
            'RLYRIKPEALKFTSIVDTPEVIQAKINAVQISEPLYRDAWEKEKANVNVPADTPLMLQSKINALQISNKRYQQAWEDVKM'
            'TGYDLRADAIGIQHAKASRDIASDYLYKTAYEKQKGHYIGCRSAKEDPKLVWAANVLKMQNDRLYKKAYNDHKAKISIPV'
            'DMVSISAAKEGQALASDVDYRHYLHRWSCFPDQNDVIQARKAYDLQSDALYKADLEWLRGIGWMPQGSPEVLRVKNAQNI'
            'FCDSVYRTPVVNLKYTSIVDTPEVVLAKSNAENISIPKYREVWDKDKTSIHIMPDTPEINLARANALNVSNKLYREGWDE'
            'MKAGCDVRLDAIPIQAAKASREIASDYKYKLDHEKQKGHYVGTLTARDDNKIRWALIADKLQNEREYRLDWAKWKAKIQS'
            'PVDMLSILHSKNSQALVSDMDYRNYLHQWTCMPDQNDVIQAKKAYELQSDNVYKADLEWLRGIGWMPNDSVSVNHAKHAA'
            'DIFSEKKYRTKIETLNFTPVDDRVDYVTAKQSGEILDDIKYRKDWNATKSKYTLTETPLLHTAQEAARILDQYLYKEGWE'
            'RQKATGYILPPDAVPFVHAHHCNDVQSELKYKAEHVKQKGHYVGVPTMRDDPKLVWFEHAGQIQNERLYKEDYHKTKAKI'
            'NIPADMVSVLAAKQGQTLVSDIDYRNYLHQWMCHPDQNDVIQARKAYDLQSDNVYRADLEWLRGIGWIPLDSVDHVRVTK'
            'NQEMMSQIKYKKNALENYPNFRSVVDPPEIVLAKINSVNQSDVKYKETFNKAKGKYTFSPDTPHISHSKDMGKLYSTILY'
            'KGAWEGTKAYGYTLDERYIPIVGAKHADLVNSELKYKETYEKQKGHYLAGKVIGEFPGVVHCLDFQKMRSALNYRKHYED'
            'TKANVHIPNDMMNHVLAKRCQYILSDLEYRHYFHQWTSLLEEPNVIRVRNAQEILSDNVYKDDLNWLKGIGCYVWDTPQI'
            'LHAKKSYDLQSQLQYTAAGKENLQNYNLVTDTPLYVTAVQSGINASEVKYKENYHQIKDKYTTVLETVDYDRTRNLKNLY'
            'SSNLYKEAWDRVKATSYILPSSTLSLTHAKNQKHLASHIKYREEYEKFKALYTLPRSVDDDPNTARCLRVGKLNIDRLYR'
            'SVYEKNKMKIHIVPDMVEMVTAKDSQKKVSEIDYRLRLHEWICHPDLQVNDHVRKVTDQISDIVYKDDLNWLKGIGCYVW'
            'DTPEILHAKHAYDLRDDIKYKAHMLKTRNDYKLVTDTPVYVQAVKSGKQLSDAVYHYDYVHSVRGKVAPTTKTVDLDRAL'
            'HAYKLQSSNLYKTSLRTLPTGYRLPGDTPHFKHIKDTRYMSSYFKYKEAYEHTKAYGYTLGPKDVPFVHVRRVNNVTSER'
            'LYRELYHKLKDKIHTTPDTPEIRQVKKTQEAVSELIYKSDFFKMQGHMISLPYTPQVIHCRYVGDITSDIKYKEDLQVLK'
            'GFGCFLYDTPDMVRSRHLRKLWSNYLYTDKARKMRDKYKVVLDTPEYRKVQELKTHLSELVYRAAGKKQKSIFTSVPDTP'
            'DLLRAKRGQKLQSQYLYVELATKERPHHHAGNQTTALKHAKDVKDMVSEKKYKIQYEKMKDKYTPVPDTPILIRAKRAYW'
            'NASDLRYKETFQKTKGKYHTVKDALDIVYHRKVTDDISKIKYKENYMSQLGIWRSIPDRPEHFHHRAVTDTVSDVKYKED'
            'LTWLKGIGCYAYDTPDFTLAEKNKTLYSKYKYKEVFERTKSDFKYVADSPINRHFKYATQLMNERKYKSSAKMFLQHGCN'
            'EILRPDMLTALYNSHMWSQIKYRKNYEKSKDKFTSIVDTPEHLRTTKVNKQISDILYKLEYNKAKPRGYTTIHDTPMLLH'
            'VRKVKDEVSDLKYKEVYQRNKSNCTIEPDAVHIKAAKDAYKVNTNLDYKKQYEANKAHWKWTPDRPDFLQAAKSSLQQSD'
            'FEYKLDREFLKGCKLSVTDDKNTVLALRNTLIESDLKYKEKHVKERGTCHAVPDTPQILLAKTVSNLVSENKYKDHVKKH'
            'LAQGSYTTLPETRDTVHVKEVTKHVSDTNYKKKFVKEKGKSNYSIMLEPPEVKHAMEVAKKQSDVAYRKDAKENLHYTTV'
            'ADRPDIKKATQAAKQASEVEYRAKHRKEGSHGLSMLGRPDIEMAKKAAKLSSQVKYRENFDKEKGKTPKYNPKDSQLYKV'
            'MKDANNLASEVKYKADLKKLHKPVTDMKESLIMNHVLNTSQLASSYQYKKKYEKSKGHYHTIPDNLEQLHLKEATELQSI'
            'VKYKEKYEKERGKPMLDFETPTYITAKESQQMQSGKEYRKDYEESIKGRNLTGLEVTPALLHVKYATKIASEKEYRKDLE'
            'ESIRGKGLTEMEDTPDMLRAKNATQILNEKEYKRDLELEVKGRGLNAMANETPDFMRARNATDIASQIKYKQSAEMEKAN'
            'FTSVVDTPEIIHAQQVKNLSSQKKYKEDAEKSMSYYETVLDTPEIQRVRENQKNFSLLQYQCDLKNSKGKITVVQDTPEI'
            'LRVKENQKNFSSVLYKEDVSPGTAIGKTPEMMRVKQTQDHISSVKYKEAIGQGTPIPDLPEVKRVKETQKHISSVMYKEN'
            'LGTGIPTTVTPEIERVKRNQENFSSVLYKENLGKGIPTPITPEMERVKRNQENFSSILYKENLSKGTPLPVTPEMERVKL'
            'NQENFSSVLYKENVGKGIPIPITPEMERVKHNQENFSSVLYKENLGTGIPIPITPEMQRVKHNQENLSSVLYKENMGKGT'
            'PLPVTPEMERVKHNQENISSVLYKENMGKGTPLPVTPEMERVKHNQENISSVLYKENMGKGTPLAVTPEMERVKHNQENI'
            'SSVLYKENVGKATATPVTPEMQRVKRNQENISSVLYKENLGKATPTPFTPEMERVKRNQENFSSVLYKENMRKATPTPVT'
            'PEMERAKRNQENISSVLYSDSFRKQIQGKAAYVLDTPEMRRVRETQRHISTVKYHEDFEKHKGCFTPVVTDPITERVKKN'
            'MQDFSDINYRGIQRKVVEMEQKRNDQDQETITGLRVWRTNPGSVFDYDPAEDNIQSRSLHMINVQAQRRSREQSRSASAL'
            'SISGGEEKSEHSEAPDHHLSTYSDGGVFAVSTAYKHAKTTELPQQRSSSVATQQTTVSSIPSHPSTAGKIFRAMYDYMAA'
            'DADEVSFKDGDAIINVQAIDEGWMYGTVQRTGRTGMLPANYVEAI*'
        )

        region_seqs = [
            'TPYIAHSQKMQDLFSPNKYKEKFEKTKG',
            'DTPELRRIKKVQDQLSEVKYRMDGD',
            'DIEHAKKVSQQVSKVLYKQNWEDTK',
            'DAPELVQAVKNTAMFSKKLYTEDWEADK',
            'DPPDIEFAKKVTNQVSKQKYKEDYEN',
            'ETPEYKMNKKAGVAASKVKYKEDYEKNKG',
            'NPQLRQLKAAGDALSDKLYKENYEKTKA',
            'DSPVLLQAQVNSKQLSDLNYKAKHESEK',
            'DTPAFIQHKVNAYNLSDNLYKQDWEKSKA',
            'DAIPLLAAKANTKNTSDVMYKKDYEKNKG',
            'DDPKMLHSLKVAKNQSDRLYKENYEKTKA',
            'DAPQFIQHRVNAYNLSDNVYKQDWEKSKA',
            'DAIPLLAAKANTKNTSDVMYKKDYEKSKG',
            'DDPKMLHSLKTAKNQSDREYRKDYEKSK',
            'DSINVDLAKKAYALQSDVEYKADYNSW',
            'DAIPIRAAKAARQAASDVQYKKDYEKAKG',
            'DDPKLVHYMNVAKIQSDREYKKDYEKTKS',
            'DAMDLELSKNMMQIQSDNVYKEDYNNWM',
            'DSPVMVQAKQNTKQVSDILYKAKGEDVKH',
            'DLPQFLQAKCNAYNISDVCYKRDWYD',
            'DAIPITAAKASRNIASDYKYKEAYEKSKG',
            'DDPKLVHYMNVAKLQSDREYKKNYENTK',
            'PQFIQAKVNALNMSDAHYKADWKKTI',
            'DAIPIVAAKSSRNIASDCKYKEAYEKAKG',
            'DDPKLVHYMNVAKIQSDREYKKGYEASK',
            'DTPDILLSRVNQITMSDKLYKAGWEEEKK',
            'DAIAIKAARASRDIASDYKYKKAYEQAKG',
            'DDPKLVHFMQVAKMQSDREYKKGYEKSK',
            'DAMSFELAKNMMQIQSDNQYKADYA',
            'DSMNMVLAQNNAKIMNEHLYKQAWEADK',
            'DIPQIILAKANAINMSDKLYKLSLEESKK',
            'DAIPIKAAKASRDIASDYKYKYNYEKGKG',
            'DDPKLVHSMQVAKMQSDREYKKNYENTK',
            'DSMDMVLAKQNAHTMNKHLYTIDWNKDK',
            'DTPDILQAKQNQTLYSQKLYKLGWEEA',
            'DAISVQLAKASRDIASDYKYKQGYRKQLG',
            'DDPKLVLSMNVAKMQSEREYKKDFEKWK',
            'QNDVVQAKKVYELQSENLYKSDLEWLRG',
            'DAMDIVLAKTNAKNRSDRLYREAWDKDK',
            'DTPDIVLAKANLINTSDKLYRMGYEELK',
            'DAIPIKAAKASREIASEYKYKEGFRKQLG',
            'DDPKMMWSMHVAKIQSDREYKKDFEKWK',
            'DVIHARQAYDLQSDNLYKSDLQWLKG',
            'DSIPMVLAKNNAITMNHRLYTEAWDKDK',
            'DTPEVLLAKQNKVNYSEKLYKLGLEEAK',
            'DDPKMMWSMHVAKIQSDREYKKDFEKWK',
            'DSLEQVLAKNNAITMNKRLYTEAWDKDK',
            'DTPEIMLARMNKINYSESLYKLANEEAK',
            'DDPKMMWSMHVAKIQSDREYKKDFEKWK',
            'DSLEQVLAKNNALNMNKRLYTEAWDKDK',
            'DTPEIMLARQNKINYSETLYKLANEEAK',
            'DAIPIVAAKASRDVISDYKYKDGYRKQLG',
            'DDPKMMWSMHVAKIQSDREYKKDFEKWK',
            'DSLEQVLAKNNALNMNKRLYTEAWDKDK',
            'DTPEIMLARQNKINYSESLYRQAMEEAKK',
            'DAIPIVAAKASRDIASDYKYKEAYRKQLG',
            'DDPKIMWSLHIAKVQSDREYKKDFEKYK',
            'QNDIIHARKAYDLQSDNLYKSDLEWMKG',
            'DTPEQVLAKNNALNMNKRLYTEAWDNDKK',
            'DTPEIMLAKLNRINYSDKLYKLALEESKK',
            'DAIPIQAAKASRDIASDYKYKEGYRKQLG',
            'DDPKMMWSIHVAKIQSDREYKKEFEKWK',
            'DTPEIVLAKNNALTMSKHLYTEAWDADK',
            'DTPDILLAKSNSANISQKLYTKGWDESK',
            'DAISIKSAKASRDIASDYKYKEAYEKQKG',
            'DDPKIMCAIHAGKIQSEREYKKEFQKWK',
            'QNDIIQAKKAYDLQSDSVYKADLEWLRG',
            'DTPEVVLAKANSLQISEKLYQEAWNKDKS',
            'DTPEMLQAHINALQISNKLYQKDWNDAK',
            'DAIEIKHAKASREIASEYKYKEGYRKQLG',
            'DDPKSVWAIHAAKIQSDREYKKAYEKSKG',
            'NDVIQAKKAYDLQSDNLYKSDLEWLKG',
            'DTPEVIQAKINAVQISEPLYRDAWEKEKA',
            'DAIGIQHAKASRDIASDYLYKTAYEKQKG',
            'NDVIQARKAYDLQSDSVYKADLEWLRG',
            'DTPEVVLAKANSLQISEKLYQEAWNKDKS',
            'DTPEMLQAHINALQISNKLYQKDWNDTK',
            'DAIEIKHAKASREIASEYKYKEGYRKQLG',
            'DDPKSVWAIHAAKIQSDREYKKAYEKSKG',
            'NDVIQAKKAYDLQSDNLYKSDLEWLKG',
            'DTPEVIQAKINAVQISEPLYRNAWEKEKA',
            'DAIGIQHAKASRDIASDYLYKTAYEKQKG',
            'NDVIQARKAYDLQSDSVYKADLEWLRG',
            'DTPEVVLAKANSLQISEKLYQEAWNKDKS',
            'DTPEMLQAHINALQISNKLYQKDWNDTK',
            'DAIEIKHAKASREIASEYKYKEGYRKQLG',
            'DDPKSVWAIHAAKIQSDREYKKAYEKSKG',
            'NDVIQAKKAYDLQSDNLYKSDLEWLKG',
            'DTPEVIQAKINAVQISEPLYRDAWEKEKA',
            'DAIGIQHAKASRDIASDYLYKTAYEKQKG',
            'NDVIQARKAYDLQSDALYKADLEWLRG',
            'DTPEVVLAKSNAENISIPKYREVWDKDK',
            'DTPEINLARANALNVSNKLYREGWDEMKA',
            'DAIPIQAAKASREIASDYKYKLDHEKQKG',
            'NDVIQAKKAYELQSDNVYKADLEWLRG',
            'DRVDYVTAKQSGEILDDIKYRKDWNATKS',
            'DAVPFVHAHHCNDVQSELKYKAEHVKQKG',
            'DDPKLVWFEHAGQIQNERLYKEDYHKTKA',
            'NDVIQARKAYDLQSDNVYRADLEWLRG',
            'DPPEIVLAKINSVNQSDVKYKETFNKAKG',
            'DTPHISHSKDMGKLYSTILYKGAWEGTKA',
            'PIVGAKHADLVNSELKYKETYEKQKG',
            'PNVIRVRNAQEILSDNVYKDDLNWLKG',
            'DTPQILHAKKSYDLQSQLQYTAAGKEN',
            'DTPLYVTAVQSGINASEVKYKENYHQIK',
            'LSLTHAKNQKHLASHIKYREEYEKFKA',
            'DTPEILHAKHAYDLRDDIKYKAH',
            'DTPHFKHIKDTRYMSSYFKYKEAYEHTKA',
            'DTPEIRQVKKTQEAVSELIYKSDFFKMQG',
            'TPQVIHCRYVGDITSDIKYKEDLQVLK',
            'DTPDMVRSRHLRKLWSNYLYTDKARKMR',
            'DTPILIRAKRAYWNASDLRYKETFQKTKG',
            'DRPEHFHHRAVTDTVSDVKYKEDLTWLKG',
            'DTPDFTLAEKNKTLYSKYKYKEVFERTKS',
            'RPDMLTALYNSHMWSQIKYRKNYEKSK',
            'DTPEHLRTTKVNKQISDILYKLEYNKAK',
            'DTPMLLHVRKVKDEVSDLKYKEVYQRNK',
            'DTPQILLAKTVSNLVSENKYKDHVKK',
            'ETRDTVHVKEVTKHVSDTNYKKKFVKEKG',
            'RPDIEMAKKAAKLSSQVKYRENFDKEKG',
            'DNLEQLHLKEATELQSIVKYKEKYEKERG',
            'ETPTYITAKESQQMQSGKEYRKDYEESI',
            'TPALLHVKYATKIASEKEYRKDLEES',
            'DTPDMLRAKNATQILNEKEYKRDLE',
            'ETPDFMRARNATDIASQIKYKQSAEMEKA',
            'DTPEIIHAQQVKNLSSQKKYKEDAEKSM',
            'DTPEIQRVRENQKNFSLLQYQCDLKNSKG',
            'DTPEILRVKENQKNFSSVLYKED',
            'TPEMMRVKQTQDHISSVKYKEA',
            'TPEIERVKRNQENFSSVLYKENLGK',
            'TPEMERVKRNQENFSSILYKENL',
            'TPEMERVKLNQENFSSVLYKEN',
            'TPEMERVKHNQENFSSVLYKEN',
            'TPEMQRVKHNQENLSSVLYKENM',
            'TPEMERVKHNQENISSVLYKENM',
            'TPEMERVKHNQENISSVLYKENM',
            'TPEMERVKHNQENISSVLYKENVGK',
            'TPEMQRVKRNQENISSVLYKENLGKA',
            'TPEMERVKRNQENFSSVLYKENMRKA',
            'TPEMERAKRNQENISSVLYSDSFRKQI',
            'DTPEMRRVRETQRHISTVKYHEDFEKHKG',
        ]
        regions = []
        p = 1
        for seq in region_seqs:
            regions.append(DomainRegion(p, p + len(seq) - 1, seq=seq))
            p += len(seq)
        d = Domain('name', regions=regions)
        with pytest.raises(UserWarning):
            d.align_seq(input_seq)


class TestCalculateORF:
    @timeout_decorator.timeout(20)
    def test_very_long(self):
        # load the sequence
        with open(get_data('calc_orf_test_sequence.fa'), 'r') as fh:
            seq = fh.readlines()[0].strip()
        calculate_orf(seq, 300)


class TestReferenceName:
    def test_naked_vs_naked_str(self):
        assert ReferenceName('1') == '1'
        assert ReferenceName('1') != '2'
        assert ReferenceName('1') == '1'
        assert ReferenceName('1') != '2'

    def test_naked_vs_prefixed_str(self):
        assert ReferenceName('1') == 'chr1'
        assert ReferenceName('1') != 'chr2'
        assert ReferenceName('1') == 'chr1'
        assert ReferenceName('1') != 'chr2'

    def test_prefixed_vs_prefixed_str(self):
        assert ReferenceName('chr1') == 'chr1'
        assert ReferenceName('chr1') != 'chr2'
        assert ReferenceName('chr1') == 'chr1'
        assert ReferenceName('chr1') != 'chr2'

    def test_prefixed_vs_naked_str(self):
        assert ReferenceName('chr1') == '1'
        assert ReferenceName('chr1') != '2'
        assert ReferenceName('chr1') == '1'

    def test_obj_comparison(self):
        r = ReferenceName('1')
        rprefix = ReferenceName('chr1')
        r2 = ReferenceName('2')
        r2prefix = ReferenceName('chr2')
        assert rprefix == r
        assert r == rprefix
        assert ReferenceName('chr1') == rprefix
        assert ReferenceName('1') == r
        assert rprefix != r2
        assert rprefix != r2prefix
        assert r != r2
        assert r != r2prefix
        assert r == rprefix
        assert r != r2prefix
        assert not r != rprefix

    def test_lt(self):
        r = ReferenceName('1')
        rprefix = ReferenceName('chr1')
        r2 = ReferenceName('2')
        r2prefix = ReferenceName('chr2')
        assert r <= rprefix
        assert not r < rprefix
        assert not rprefix < r
        assert rprefix <= r
        for chr1, chr2 in itertools.product([r, rprefix], [r2, r2prefix]):
            assert chr1 < chr2
            assert chr1 <= chr2

    def test_alpha_sort(self):
        assert ReferenceName('10') < ReferenceName('3')
        assert ReferenceName('10') < ReferenceName('chr3')
        assert ReferenceName('chr10') < ReferenceName('3')
        assert ReferenceName('chr10') < ReferenceName('chr3')

    def test_gt(self):
        r = ReferenceName('1')
        rprefix = ReferenceName('chr1')
        r2 = ReferenceName('2')
        r2prefix = ReferenceName('chr2')
        assert rprefix >= r
        assert r >= rprefix
        assert not r > rprefix
        assert not rprefix > r
        for chr1, chr2 in itertools.product([r, rprefix], [r2, r2prefix]):
            assert chr2 > chr1
            assert chr2 >= chr1

    def test_hash(self):
        assert ReferenceName('3') in {ReferenceName('3')}
        assert ReferenceName('3') in {ReferenceName('chr3')}


class TestIndelCall:
    def test_duplication_in_repeat(self):
        ref = 'ASFHGHGSFSFSLLLLLL' 'FLLLLSFSLMVPWSFKW'
        mut = 'ASFHGHGSFSFSLLLLLLL' 'FLLLLSFSLMVPWSFKW'

        call = IndelCall(ref, mut)
        print(call)

        assert call.nterm_aligned == 18
        assert call.cterm_aligned == len(ref) - 13 + 1
        assert call.is_dup

        assert call.hgvs_protein_notation() == 'p.L18dupL'

    def test_nterminal_extension(self):

        ref = 'MABCDEFGH'
        mut = 'MAFMABCDEFGH'

        call = IndelCall(ref, mut)
        print(call)
        assert not call.nterm_aligned
        assert call.cterm_aligned == len(call.ref_seq) - 1 + 1
        assert not call.is_dup
        assert call.ins_seq == 'MAF'
        assert call.del_seq == ''

        assert call.hgvs_protein_notation() == 'p.M1ext-3'

    def test_nterminal_deletion(self):
        ref = 'MABCDEFGH'
        mut = 'CDEFGH'

        call = IndelCall(ref, mut)
        print(call)
        assert not call.nterm_aligned
        assert call.cterm_aligned == len(call.ref_seq) - 4 + 1
        assert not call.is_dup
        assert call.ins_seq == ''
        assert call.del_seq == 'MAB'

        assert call.hgvs_protein_notation() == 'p.M1_B3delMAB'

    def test_cterminal_deletion(self):
        ref = 'MABCDEFGH'
        mut = 'MABCDE'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 6
        assert not call.cterm_aligned
        assert not call.is_dup
        assert call.ins_seq == ''
        assert call.del_seq == 'FGH'

        assert call.hgvs_protein_notation() == 'p.F7_H9delFGH'

    def test_cterminal_extension(self):
        ref = 'MABCDEFGH'
        mut = 'MABCDEFGHIJK'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 9
        assert not call.cterm_aligned
        assert not call.is_dup
        assert call.ins_seq == 'IJK'
        assert call.del_seq == ''

        assert call.hgvs_protein_notation() == 'p.H9ext3'

    def test_cterminal_stop_extension(self):
        ref = 'MABCDEFGH*'
        mut = 'MABCDEFGHIJK*'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 9
        assert not call.cterm_aligned
        assert not call.is_dup
        assert call.ins_seq == 'IJK'
        assert call.del_seq == ''

        assert call.hgvs_protein_notation() == 'p.*10ext*3'

    def test_cterminal_no_orf_ext(self):
        ref = 'MABCDEFGH'
        mut = 'MABCDEFGHIJK*'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 9
        assert not call.cterm_aligned
        assert not call.is_dup
        assert call.ins_seq == 'IJK*'
        assert call.del_seq == ''

        assert call.hgvs_protein_notation() == 'p.H9ext*4'

    def test_single_aa_insertion(self):
        ref = 'MABCDEFGH'
        mut = 'MABCKDEFGH'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 4
        assert call.cterm_aligned == len(call.ref_seq) - 5 + 1
        assert not call.is_dup
        assert call.ins_seq == 'K'
        assert call.del_seq == ''

        assert call.hgvs_protein_notation() == 'p.C4_D5insK'

    def test_insertion(self):
        ref = 'MABCDEFGH'
        mut = 'MABCKADEFGH'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 4
        assert call.cterm_aligned == len(call.ref_seq) - 5 + 1
        assert not call.is_dup
        assert call.ins_seq == 'KA'
        assert call.del_seq == ''

        assert call.hgvs_protein_notation() == 'p.C4_D5insKA'

    def test_single_aa_deletion(self):
        ref = 'MABCDEFGH'
        mut = 'MABCEFGH'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 4
        assert call.cterm_aligned == len(call.ref_seq) - 6 + 1
        assert not call.is_dup
        assert call.ins_seq == ''
        assert call.del_seq == 'D'

        assert call.hgvs_protein_notation() == 'p.D5delD'

    def test_deletion(self):
        ref = 'MABCDEFGH'
        mut = 'MABCFGH'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 4
        assert call.cterm_aligned == len(call.ref_seq) - 7 + 1
        assert not call.is_dup
        assert call.ins_seq == ''
        assert call.del_seq == 'DE'

        assert call.hgvs_protein_notation() == 'p.D5_E6delDE'

    def test_deletion_in_repeat(self):
        ref = 'MABCDEEEEEEFGH'
        mut = 'MABCDEEEEFGH'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 9
        assert call.cterm_aligned == len(call.ref_seq) - 8 + 1
        assert not call.is_dup
        assert call.ins_seq == ''
        assert call.del_seq == 'EE'

        assert call.hgvs_protein_notation() == 'p.E10_E11delEE'

    def test_insertion_in_repeat(self):
        ref = 'MABCDEEEEFGH'
        mut = 'MABCDEEEEEEFGH'

        call = IndelCall(ref, mut)
        print(call)
        assert call.nterm_aligned == 9
        assert call.cterm_aligned == len(call.ref_seq) - 6 + 1
        assert call.is_dup
        assert call.ins_seq == 'EE'
        assert call.del_seq == ''

        assert call.hgvs_protein_notation() == 'p.E8_E9dupEE'
