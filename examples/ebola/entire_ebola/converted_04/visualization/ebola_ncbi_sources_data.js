  var sources = { data: [
{ id: 0,	name:'eboVir3',	title: 'G3686v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 1,	name:'KJ660346v2',	title: 'Guinea_Kissidougou-C15_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 2,	name:'KJ660347v2',	title: 'Guinea_Gueckedou-C07_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 3,	name:'KJ660348v2',	title: 'Guinea_Gueckedou-C05_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 4,	name:'KM233050v1',	title: 'G3713v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 5,	name:'KM233069v1',	title: 'G3770v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 6,	name:'KM233070v1',	title: 'G3770v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 7,	name:'KM233104v1',	title: 'G3834_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18921}
,{ id: 8,	name:'KM233109v1',	title: 'G3846_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 9,	name:'KM233113v1',	title: 'G3856v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 10,	name:'AF086833v2',	title: 'AF086833v2_1976',	 group_name: 'Zaire(DRC) 1976-7',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 11,	name:'AF272001v1',	title: 'Mayinga_1976',	 group_name: 'Zaire(DRC) 1976-7',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 12,	name:'AY142960v1',	title: 'Mayinga_2002',	 group_name: 'Zaire(DRC) 1976-7',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 13,	name:'EU224440v2',	title: 'GuineaPig_Mayinga_2007',	 group_name: 'Zaire(DRC) 1976-7',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 14,	name:'KC242791v1',	title: 'Bonduni_1977',	 group_name: 'Zaire(DRC) 1976-7',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 15,	name:'KC242792v1',	title: 'Gabon_1994',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 16,	name:'KC242794v1',	title: '2Nza_1996',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 17,	name:'KC242796v1',	title: '13625Kikwit_1995',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 18,	name:'KC242798v1',	title: '1Ikot_Gabon_1996',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 19,	name:'KC242799v1',	title: '13709Kikwit_1995',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 20,	name:'KC242801v1',	title: 'deRoover_1976',	 group_name: 'Zaire(DRC) 1976-7',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 21,	name:'NC_002549v1',	title: 'NC_002549v1_1976',	 group_name: 'Zaire(DRC) 1976-7',	bundle_ID: 0,	weight: -1,	length: 18957}
,{ id: 22,	name:'AY354458v1',	title: 'Zaire_1995',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18958}
,{ id: 23,	name:'KC242784v1',	title: 'Luebo9_2007',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 24,	name:'KC242785v1',	title: 'Luebo0_2007',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 25,	name:'KC242786v1',	title: 'Luebo1_2007',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 26,	name:'KC242787v1',	title: 'Luebo23_2007',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 27,	name:'KC242788v1',	title: 'Luebo43_2007',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 28,	name:'KC242789v1',	title: 'Luebo4_2007',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 29,	name:'KC242790v1',	title: 'Luebo5_2007',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 30,	name:'KC242793v1',	title: '1Eko_1996',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 31,	name:'KC242795v1',	title: '1Mbie_Gabon_1996',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 32,	name:'KC242797v1',	title: '1Oba_Gabon_1996',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 33,	name:'KC242800v1',	title: 'Ilembe_2002',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 34,	name:'AF499101v1',	title: 'Mouse_Mayinga_2002',	 group_name: 'Zaire(DRC) 1976-7',	bundle_ID: 0,	weight: -1,	length: 18958}
,{ id: 35,	name:'FJ217162v1',	title: 'Cote_dIvoire_CIEBOV_1994',	 group_name: 'Bundibugyo 2007',	bundle_ID: -1,	weight: -1,	length: 17983}
,{ id: 36,	name:'NC_014372v1',	title: 'Cote_dIvoire_1994',	 group_name: 'Bundibugyo 2007',	bundle_ID: -1,	weight: -1,	length: 17983}
,{ id: 37,	name:'FJ217161v1',	title: 'Bundibugyo_Uganda_2007',	 group_name: 'Bundibugyo 2007',	bundle_ID: 3,	weight: -1,	length: 18606}
,{ id: 38,	name:'NC_014373v1',	title: 'Bundibugyo_2007',	 group_name: 'Bundibugyo 2007',	bundle_ID: 3,	weight: -1,	length: 18606}
,{ id: 39,	name:'KC545395v1',	title: 'EboBund-122_2012',	 group_name: 'Bundibugyo 2007',	bundle_ID: 3,	weight: 0,	length: 18484}
,{ id: 40,	name:'KC545394v1',	title: 'EboBund-120_2012',	 group_name: 'Bundibugyo 2007',	bundle_ID: 3,	weight: 50,	length: 18489}
,{ id: 41,	name:'KC545393v1',	title: 'EboBund-112_2012',	 group_name: 'Bundibugyo 2007',	bundle_ID: 3,	weight: 100,	length: 18479}
,{ id: 42,	name:'KC545396v1',	title: 'EboBund-14_2012',	 group_name: 'Bundibugyo 2007',	bundle_ID: 3,	weight: 99,	length: 18479}
,{ id: 43,	name:'FJ621584v1',	title: 'Reston08-C_2008',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 17684}
,{ id: 44,	name:'JX477166v1',	title: 'Alice_TX_USA_MkCQ8167_1996',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 18570}
,{ id: 45,	name:'AY769362v1',	title: 'reconstructReston_2008',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 18489}
,{ id: 46,	name:'EU338380v1',	title: 'Yambio_2004',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18736}
,{ id: 47,	name:'KC242783v2',	title: 'Maleo_1979',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18738}
,{ id: 48,	name:'JX477165v1',	title: 'Reston09-A_2009',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 18614}
,{ id: 49,	name:'AF522874v1',	title: 'Reston_PA_1990',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 18461}
,{ id: 50,	name:'NC_004161v1',	title: 'Pennsylvania_1990',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 18461}
,{ id: 51,	name:'FJ621583v1',	title: 'Reston08-A_2008',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 18548}
,{ id: 52,	name:'KC589025v1',	title: 'EboSud-639_2012',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18468}
,{ id: 53,	name:'FJ968794v1',	title: 'Boniface_1976',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18731}
,{ id: 54,	name:'AY729654v1',	title: 'Gulu_Uganda_2000',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18591}
,{ id: 55,	name:'NC_006432v1',	title: 'Gulu_2000',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18592}
,{ id: 56,	name:'KC545389v1',	title: 'EboSud-602_2012',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18650}
,{ id: 57,	name:'KC545390v1',	title: 'EboSud-603_2012',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18650}
,{ id: 58,	name:'KC545391v1',	title: 'EboSud-609_2012',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18650}
,{ id: 59,	name:'KC545392v1',	title: 'EboSud-682_2012',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18650}
,{ id: 60,	name:'JN638998v1',	title: 'Nakisamata_2011',	 group_name: 'Sudan 1976-9',	bundle_ID: 1,	weight: -1,	length: 18587}
,{ id: 61,	name:'KM233053v1',	title: 'G3724_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 62,	name:'KM233056v1',	title: 'G3735v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18956}
,{ id: 63,	name:'AB050936v1',	title: 'Reston_1996',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 18555}
,{ id: 64,	name:'KM034557v1',	title: 'G3677v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18955}
,{ id: 65,	name:'KM233045v1',	title: 'EM124v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18919}
,{ id: 66,	name:'KM233110v1',	title: 'G3848_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18955}
,{ id: 67,	name:'KM233096v1',	title: 'G3822_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18953}
,{ id: 68,	name:'KM233057v1',	title: 'G3735v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18953}
,{ id: 69,	name:'KM233063v1',	title: 'G3764_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18953}
,{ id: 70,	name:'KM233089v1',	title: 'G3814_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18914}
,{ id: 71,	name:'KM233051v1',	title: 'G3713v3_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18951}
,{ id: 72,	name:'KM233097v1',	title: 'G3823_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18952}
,{ id: 73,	name:'KM034560v1',	title: 'G3682v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18951}
,{ id: 74,	name:'KM233039v1',	title: 'EM112_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18951}
,{ id: 75,	name:'KM233043v1',	title: 'EM120_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18915}
,{ id: 76,	name:'KM233099v1',	title: 'G3825v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18951}
,{ id: 77,	name:'KM034554v1',	title: 'G3676v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18948}
,{ id: 78,	name:'KM034555v1',	title: 'G3676v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18949}
,{ id: 79,	name:'KM233072v1',	title: 'G3782_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18949}
,{ id: 80,	name:'KM233092v1',	title: 'G3818_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18948}
,{ id: 81,	name:'KM233098v1',	title: 'G3825v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18949}
,{ id: 82,	name:'KM233103v1',	title: 'G3831_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18949}
,{ id: 83,	name:'KM233049v1',	title: 'G3707_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18946}
,{ id: 84,	name:'KM233054v1',	title: 'G3729_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18946}
,{ id: 85,	name:'KM233080v1',	title: 'G3799_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18945}
,{ id: 86,	name:'KM233102v1',	title: 'G3829_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18944}
,{ id: 87,	name:'KM233106v1',	title: 'G3840_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18944}
,{ id: 88,	name:'JQ352763v1',	title: 'Kikwit_1995',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18940}
,{ id: 89,	name:'KM233052v1',	title: 'G3713v4_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18941}
,{ id: 90,	name:'KM233114v1',	title: 'G3856v3_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18939}
,{ id: 91,	name:'KM233062v1',	title: 'G3758_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18937}
,{ id: 92,	name:'KM034556v1',	title: 'G3677v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18932}
,{ id: 93,	name:'KM233077v1',	title: 'G3795_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18897}
,{ id: 94,	name:'KM233100v1',	title: 'G3826_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18932}
,{ id: 95,	name:'KM233075v1',	title: 'G3788_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18931}
,{ id: 96,	name:'KM034553v1',	title: 'G3670v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18927}
,{ id: 97,	name:'KM034558v1',	title: 'G3679v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18915}
,{ id: 98,	name:'KM034551v1',	title: 'EM096_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18926}
,{ id: 99,	name:'KM034559v1',	title: 'G3680v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18926}
,{ id: 100,	name:'KM233093v1',	title: 'G3819_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18880}
,{ id: 101,	name:'KM233116v1',	title: 'NM042v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18926}
,{ id: 102,	name:'KM233061v1',	title: 'G3752_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18890}
,{ id: 103,	name:'KM233065v1',	title: 'G3769v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18925}
,{ id: 104,	name:'KM233068v1',	title: 'G3769v4_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18925}
,{ id: 105,	name:'KM233074v1',	title: 'G3787_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18890}
,{ id: 106,	name:'KM233081v1',	title: 'G3800_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18925}
,{ id: 107,	name:'KM233094v1',	title: 'G3820_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18925}
,{ id: 108,	name:'KM034561v1',	title: 'G3683v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18908}
,{ id: 109,	name:'KM233040v1',	title: 'EM113_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18922}
,{ id: 110,	name:'KM233041v1',	title: 'EM115_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18886}
,{ id: 111,	name:'KM233047v1',	title: 'EM124v3_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18921}
,{ id: 112,	name:'KM233084v1',	title: 'G3807_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18921}
,{ id: 113,	name:'KM233088v1',	title: 'G3810v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18921}
,{ id: 114,	name:'KM233112v1',	title: 'G3851_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18922}
,{ id: 115,	name:'KM233091v1',	title: 'G3817_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18886}
,{ id: 116,	name:'KM233085v1',	title: 'G3808_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18882}
,{ id: 117,	name:'FJ621585v1',	title: 'Reston08-E_2008',	 group_name: 'Reston 1989-90',	bundle_ID: 2,	weight: -1,	length: 17678}
,{ id: 118,	name:'KM233046v1',	title: 'EM124v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18867}
,{ id: 119,	name:'KM233055v1',	title: 'G3734v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18878}
,{ id: 120,	name:'KM233035v1',	title: 'EM104_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18912}
,{ id: 121,	name:'KM233037v1',	title: 'EM110_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18913}
,{ id: 122,	name:'KM233038v1',	title: 'EM111_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18913}
,{ id: 123,	name:'KM233042v1',	title: 'EM119_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18912}
,{ id: 124,	name:'KM233044v1',	title: 'EM121_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18913}
,{ id: 125,	name:'KM233073v1',	title: 'G3786_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18912}
,{ id: 126,	name:'KM233078v1',	title: 'G3796_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18912}
,{ id: 127,	name:'KM233105v1',	title: 'G3838_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18913}
,{ id: 128,	name:'KM233108v1',	title: 'G3845_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18913}
,{ id: 129,	name:'KM233066v1',	title: 'G3769v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18876}
,{ id: 130,	name:'KM233101v1',	title: 'G3827_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18876}
,{ id: 131,	name:'KM233107v1',	title: 'G3841_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18876}
,{ id: 132,	name:'KM233115v1',	title: 'G3857_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18835}
,{ id: 133,	name:'KM034549v1',	title: 'EM095B_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18803}
,{ id: 134,	name:'KM034550v1',	title: 'EM095_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18889}
,{ id: 135,	name:'KM233067v1',	title: 'G3769v3_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18910}
,{ id: 136,	name:'KM233036v1',	title: 'EM106_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18873}
,{ id: 137,	name:'KM233064v1',	title: 'G3765v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18803}
,{ id: 138,	name:'KM233095v1',	title: 'G3821_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18862}
,{ id: 139,	name:'HQ613403v1',	title: 'M-M_2007',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18807}
,{ id: 140,	name:'KM233048v1',	title: 'EM124v4_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18829}
,{ id: 141,	name:'KM233076v1',	title: 'G3789v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18811}
,{ id: 142,	name:'KM233079v1',	title: 'G3798_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18870}
,{ id: 143,	name:'KM233086v1',	title: 'G3809_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18830}
,{ id: 144,	name:'KM233087v1',	title: 'G3810v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18859}
,{ id: 145,	name:'KM233058v1',	title: 'G3750v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18850}
,{ id: 146,	name:'KM233071v1',	title: 'G3771_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18901}
,{ id: 147,	name:'KM233059v1',	title: 'G3750v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18791}
,{ id: 148,	name:'KM233082v1',	title: 'G3805v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18812}
,{ id: 149,	name:'HQ613402v1',	title: '034-KS_2008',	 group_name: 'DRC 2007',	bundle_ID: 0,	weight: -1,	length: 18774}
,{ id: 150,	name:'KM233111v1',	title: 'G3850_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18749}
,{ id: 151,	name:'KM233090v1',	title: 'G3816_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18799}
,{ id: 152,	name:'KM233060v1',	title: 'G3750v3_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18711}
,{ id: 153,	name:'KM233083v1',	title: 'G3805v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18773}
,{ id: 154,	name:'KM233118v1',	title: 'NM042v3_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18613}
,{ id: 155,	name:'KM034552v1',	title: 'EM098_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18769}
,{ id: 156,	name:'KM233117v1',	title: 'NM042v2_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18685}
,{ id: 157,	name:'KM034563v1',	title: 'G3687v1_2014',	 group_name: 'Ebola 2014',	bundle_ID: 0,	weight: -1,	length: 18511}
,{ id: 158,	name:'NC_024781v1',	title: 'Marburg_KitumCave_Kenya_1987',	 group_name: 'Marburg 1987',	bundle_ID: -1,	weight: -1,	length: 12846}
,{ id: 159,	name:'NC_001608v3',	title: 'Marburg_MtElgon_Musoke_Kenya_1980',	 group_name: 'Marburg 1987',	bundle_ID: -1,	weight: -1,	length: 13530}]};

                    var consensuses = {data: [
{ id: 0,	name:'CONSENS0',	title: 'consensus produced by heaviest_bundle, containing 17 seqs',	sources_compatibility: [0.9994, 0.9994, 0.9993, 0.9992, 0.9996, 0.9996, 0.9996, 0.9996, 0.9996, 0.9996, 0.9702, 0.9698, 0.9702, 0.9694, 0.9702, 0.966, 0.9653, 0.9684, 0.9653, 0.9685, 0.9702, 0.9702, 0.9682, 0.9711, 0.9711, 0.9714, 0.9713, 0.9709, 0.9712, 0.9711, 0.9658, 0.9658, 0.9658, 0.9681, 0.9697, 0.682, 0.682, 0.6715, 0.6713, 0.6736, 0.6748, 0.6746, 0.6745, 0.6522, 0.6477, 0.6467, 0.6408, 0.6403, 0.6469, 0.6477, 0.6477, 0.6481, 0.6444, 0.641, 0.6453, 0.6451, 0.6454, 0.6455, 0.6454, 0.6454, 0.6425, 0.9997, 0.9997, 0.6476, 0.9998, 0.9998, 0.9997, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 1.0, 0.9999, 0.9999, 0.9999, 0.9998, 0.9998, 1.0, 0.9998, 0.9999, 0.9999, 0.9999, 1.0, 0.9999, 0.9998, 0.9999, 0.9685, 0.9999, 0.9999, 1.0, 1.0, 1.0, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9997, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 1.0, 1.0, 0.9997, 0.9999, 0.9999, 0.9999, 1.0, 0.9999, 0.9999, 0.9998, 1.0, 0.651, 0.9999, 1.0, 0.9999, 0.9999, 0.9999, 0.9998, 0.9999, 1.0, 1.0, 1.0, 0.9999, 0.9999, 0.9999, 1.0, 0.9999, 0.9998, 0.9998, 0.9999, 0.9999, 0.9999, 0.9999, 0.9715, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 0.971, 0.9999, 0.9999, 0.9999, 0.9999, 0.9999, 1.0, 0.9999, 0.9973, 0.5571, 0.5565],	length: 18952 }
,{ id: 0,	name:'CONSENS1',	title: 'consensus produced by heaviest_bundle, containing 4 seqs',	sources_compatibility: [0.6353, 0.6353, 0.6353, 0.6351, 0.6352, 0.6352, 0.6352, 0.6349, 0.6352, 0.6352, 0.6311, 0.6309, 0.6311, 0.6306, 0.6312, 0.6306, 0.6301, 0.6313, 0.6301, 0.6313, 0.6311, 0.6311, 0.6313, 0.6317, 0.6316, 0.6316, 0.6316, 0.6313, 0.6316, 0.6314, 0.6307, 0.6307, 0.6307, 0.63, 0.6308, 0.5901, 0.5901, 0.5826, 0.5825, 0.584, 0.5845, 0.5847, 0.5846, 0.5949, 0.5817, 0.5835, 0.866, 0.8651, 0.5836, 0.5856, 0.5856, 0.5849, 0.9556, 0.8815, 0.9726, 0.9728, 1.0, 0.9999, 0.9999, 0.9999, 0.9749, 0.6352, 0.6353, 0.5817, 0.6352, 0.6349, 0.6352, 0.6352, 0.6352, 0.6352, 0.6349, 0.6351, 0.6351, 0.6351, 0.6352, 0.6347, 0.6351, 0.6351, 0.6351, 0.6351, 0.6351, 0.6351, 0.6351, 0.6351, 0.635, 0.6349, 0.635, 0.635, 0.631, 0.6349, 0.6349, 0.6348, 0.6349, 0.6345, 0.6349, 0.6348, 0.6348, 0.6344, 0.6347, 0.6348, 0.6343, 0.6348, 0.6345, 0.6348, 0.6348, 0.6345, 0.6348, 0.6348, 0.6345, 0.6349, 0.6346, 0.6349, 0.6349, 0.6349, 0.6349, 0.6346, 0.6345, 0.5937, 0.6343, 0.6344, 0.6348, 0.6347, 0.6347, 0.6347, 0.6347, 0.6347, 0.6347, 0.6347, 0.6347, 0.6343, 0.6344, 0.6344, 0.6341, 0.6346, 0.635, 0.6346, 0.6343, 0.6344, 0.6342, 0.6307, 0.634, 0.6342, 0.6343, 0.6341, 0.6341, 0.6341, 0.6345, 0.6337, 0.6337, 0.6303, 0.634, 0.6339, 0.6342, 0.6338, 0.6354, 0.6338, 0.6345, 0.6344, 0.5022, 0.4916],	length: 18650 }
,{ id: 0,	name:'CONSENS2',	title: 'consensus produced by heaviest_bundle, containing 4 seqs',	sources_compatibility: [0.6312, 0.6312, 0.6312, 0.6312, 0.6311, 0.631, 0.631, 0.6308, 0.631, 0.631, 0.6304, 0.6302, 0.6304, 0.6304, 0.6304, 0.6291, 0.6287, 0.6298, 0.6287, 0.6297, 0.6303, 0.6304, 0.6298, 0.6301, 0.6301, 0.6302, 0.6302, 0.6299, 0.6303, 0.6301, 0.6291, 0.6291, 0.6291, 0.6288, 0.6301, 0.5946, 0.5946, 0.5828, 0.5828, 0.5871, 0.5877, 0.5877, 0.5875, 0.9147, 0.9314, 0.9681, 0.5787, 0.5783, 0.9173, 1.0, 1.0, 0.925, 0.5836, 0.5799, 0.5812, 0.5811, 0.5796, 0.5796, 0.5796, 0.5796, 0.5794, 0.6311, 0.631, 0.9314, 0.631, 0.6308, 0.631, 0.6309, 0.631, 0.631, 0.6307, 0.631, 0.631, 0.6309, 0.6309, 0.6307, 0.6309, 0.631, 0.631, 0.6309, 0.6308, 0.6309, 0.631, 0.6308, 0.6308, 0.6308, 0.6308, 0.6308, 0.6294, 0.6308, 0.6307, 0.6307, 0.6307, 0.6305, 0.6307, 0.6307, 0.6306, 0.6305, 0.6306, 0.6307, 0.6301, 0.6306, 0.6304, 0.6305, 0.6305, 0.6305, 0.6306, 0.6306, 0.6305, 0.6306, 0.6305, 0.6307, 0.6307, 0.6306, 0.6307, 0.6304, 0.6304, 0.9176, 0.6301, 0.6303, 0.6304, 0.6305, 0.6305, 0.6304, 0.6306, 0.6305, 0.6305, 0.6305, 0.6305, 0.6302, 0.6303, 0.6303, 0.6302, 0.6305, 0.6308, 0.6304, 0.6303, 0.6304, 0.63, 0.6296, 0.6301, 0.6303, 0.6302, 0.6301, 0.6299, 0.6299, 0.6303, 0.6297, 0.6298, 0.6293, 0.6301, 0.6296, 0.6301, 0.6297, 0.6303, 0.6297, 0.6299, 0.6304, 0.4989, 0.4877],	length: 18461 }
,{ id: 0,	name:'CONSENS3',	title: 'consensus produced by heaviest_bundle, containing 4 seqs',	sources_compatibility: [0.658, 0.6582, 0.6582, 0.6581, 0.6581, 0.658, 0.658, 0.6579, 0.658, 0.658, 0.6557, 0.6555, 0.6557, 0.6555, 0.6558, 0.6553, 0.655, 0.6562, 0.655, 0.6562, 0.6557, 0.6557, 0.656, 0.6565, 0.6564, 0.6565, 0.6565, 0.6564, 0.6564, 0.6564, 0.6554, 0.6554, 0.6554, 0.6544, 0.6555, 0.6711, 0.6711, 0.9314, 0.9316, 0.9865, 0.9934, 0.9999, 0.9997, 0.5966, 0.5862, 0.5887, 0.5769, 0.5762, 0.5846, 0.5883, 0.5883, 0.5865, 0.5831, 0.577, 0.5801, 0.5799, 0.5793, 0.5792, 0.5794, 0.5792, 0.5777, 0.6581, 0.658, 0.5857, 0.6579, 0.6579, 0.658, 0.6579, 0.6579, 0.6579, 0.6578, 0.658, 0.6579, 0.6579, 0.6579, 0.6577, 0.6579, 0.658, 0.658, 0.6578, 0.6578, 0.6579, 0.6579, 0.6578, 0.6578, 0.6578, 0.6577, 0.6578, 0.6557, 0.6578, 0.6576, 0.6577, 0.6576, 0.6575, 0.6577, 0.6575, 0.6575, 0.6574, 0.6576, 0.6577, 0.6573, 0.6576, 0.6575, 0.6576, 0.6576, 0.6575, 0.6576, 0.6576, 0.6576, 0.6577, 0.6576, 0.6577, 0.6577, 0.6577, 0.6578, 0.6576, 0.6576, 0.5952, 0.6574, 0.6575, 0.6576, 0.6576, 0.6576, 0.6576, 0.6575, 0.6576, 0.6576, 0.6576, 0.6576, 0.6575, 0.6575, 0.6574, 0.6576, 0.6579, 0.6578, 0.6576, 0.6575, 0.6578, 0.6574, 0.6565, 0.6575, 0.6576, 0.6574, 0.6575, 0.6572, 0.6573, 0.6575, 0.6572, 0.6569, 0.6558, 0.6577, 0.6574, 0.658, 0.6573, 0.6591, 0.6573, 0.6583, 0.6581, 0.5085, 0.5007],	length: 18479 }]};

                    var poagraph = { nodes_count: 46638,	sources_count:160,	mean_sources_per_node: 64.36,	consensuses_count: 4};
