from coffea.processor import Runner, FuturesExecutor
from coffea.nanoevents import NanoAODSchema
# from Custom_Processor import HiggsAnalysisProcessor
from WH_Reco_processor import HiggsAnalysisProcessor
import awkward as ak
import pyarrow.parquet as pq
import pyarrow as pa
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from coffea import processor
import glob


fileset={
  "M20_RunIISummer20UL18NanoAODv9": [
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/03AB06FF-95D3-F045-BF78-32382A914871.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/04D5460E-C7D2-E14F-A87F-06F077F6751E.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/05AFA710-C3DB-8E4F-B178-432DF57B5273.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/088028C9-1205-D548-AB0E-26486C92BD0A.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/1157DF00-1F65-2545-9D5E-4FE6F8D458E9.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/20A20C6B-1489-0443-9E01-9F5703E974B7.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/22692812-7544-D149-87EE-346B1501B9A4.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/2B3AD269-9B0A-6F46-AF95-DCF13635E4A4.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/2BC9D163-6D59-524F-9A46-7EAF7261CA07.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/354B0999-7727-9D49-A549-7FFA9800596C.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/3AE508DA-8315-4F4F-AE0C-0DDC14920F8F.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/3B6D928B-B01C-8043-A303-31FA6E9D9F78.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/3C71EA6B-BBE4-6D44-B7F1-C812611BFE63.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/4383E866-1CF6-0A4A-B231-245449C05891.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/4F784984-8AA7-B344-A940-98AD4A4D1237.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/5411CD57-C701-6243-82DD-69CF81EFAB75.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/638E0BE4-AD00-3A44-AEAD-D7FC2A9A2620.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/6638378F-E911-F34F-880B-FF5FCEBFB74B.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/762500E3-4A11-D048-9F7B-5884B11D91B8.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/97B58431-B873-AB43-AA8F-5D36EB56FEA3.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/9B889D20-431C-7046-B9F0-A72D9B3CDCDF.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/9CDD1C85-65CE-4E4A-BDDF-44D37AEFACC1.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/A1163943-3DCB-684F-BA5F-34BFAB8590BA.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/A389CCC4-CE66-7246-99C0-316D3C64DD42.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/A86A9817-C428-F642-84E4-110E7346BA26.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/A8C9D620-9344-8844-9B05-F2E27B20150C.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/AC97149D-F43B-A94B-81CD-37CCE6618A6A.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/AE8B6AFC-F5BC-E449-9926-809D8F2CFEE6.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/B2A39EB3-3950-A04F-A898-41BC0A7640EB.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/B5AE5C67-FC1A-E941-AE21-9045BF1372CB.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/B77DA2B5-9E7B-5A4B-9403-6C73CECE3CDD.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/B7D21F37-159B-4F42-B1A8-97F4732749B9.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/C2B92AC8-009F-F643-A697-4EB7D21AB0EA.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/C655E6E0-B8C2-AB49-A804-4CC8FBB45BC7.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/C94BA8CF-6A55-2A4C-8B81-6D5EEAD9BBDE.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/CBD4D93E-6303-8744-AEEF-4B65C6132F68.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/D19C3D07-9C6D-9248-8C88-A5BE008068ED.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/D1BFF410-71BA-454A-A00D-FF11F5F6BBA3.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/D27D8617-953E-304C-BF3F-C9F0D24B0B3E.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/D5566C84-8BBA-2447-BAA4-18280D90D2D8.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/D6A8A4FD-87E2-E947-9636-503D13B214A2.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/D895DE54-BC83-1D41-9E07-E2455D75EA3C.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/DCA57DD2-9AE4-D743-A489-ED76B612C2B3.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/DF50F3C2-FF86-6E44-9963-8214B66A0047.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/E188DE51-497D-424C-8D9C-D5EC5ECECFB3.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/E78DC30D-9EAF-4D4D-89A1-9F1F88153D14.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/EB0CF31A-05EB-0C48-969E-5CA4D457FEB9.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/EC1D6162-D0DF-5847-931A-398A54B4B65F.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/EDB05A41-EC07-CB4A-88B8-D77598E64852.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/F49B5445-7F6D-B442-85F9-74F5B4CF6E9D.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/F7709D83-7801-4449-97C9-888B33BA7E63.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/F7792587-54C6-F849-874D-372DE04E0D6C.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/F9035AF5-AA0A-6044-ABF1-30BC780A82E3.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/F98D73B4-4E91-4644-89DF-B48DB8519DB1.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M20-RunIISummer20UL18MiniAODv2/FD003A22-26DE-3340-8693-0C2021A54446.root"


  ],
  "M60_RunIISummer20UL18NanoAODv9": [
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/000BDF98-C383-104A-BC2E-2ACEB5C15C94.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/0311C9F9-8E1E-E647-B7ED-DB640C4FFB38.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/172C44EB-DA3F-8949-91BA-640D02971FF5.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/17ADBB6D-E775-874A-B13E-41A62D5E3999.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/201956DB-A513-DC4A-9C24-BB7977B24F52.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/34FE7756-2B2D-2F4F-8A09-EC0199CEDAAF.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/3A8D5171-432A-634B-A31B-E2149A01CFFE.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/4463B929-B03A-E540-9539-A3D02743AFCB.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/51A19845-88EA-D14A-AFBD-4C02602C9DFF.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/5EDEAD04-F4AA-1D44-8094-687E5B0B8D35.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/74985E67-9408-1649-BE9E-0981B2EE102F.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/762F7BBA-C325-C348-9449-FCD5BFD726E9.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/8413E724-9F72-7744-BA04-19C44838CD14.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/901E6E31-4CFF-A34A-96C2-EEC9D564CDF8.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/97C4EC2E-AB6B-E149-9551-198AC839196F.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/A2F42F89-6BB9-E44C-88BB-667906A42C19.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/B212E94E-42AD-F04B-9DF8-8A1654CD805B.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/BBA9736E-B02C-984F-8341-4922AE305219.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/BEF3B897-4976-4B46-A04C-F53E634DE2DC.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/C2EAA6B3-3CB9-2A44-8AD7-63B8A9FA778C.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/C5668E37-521E-CB4B-9C7F-4765172B5EAF.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/CE3BF102-2776-714E-ACA0-78B807706A83.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/D8D94DA8-B1D7-B341-90F7-789F34187864.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/D9891662-E418-CB48-9283-69A0224F35BC.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/DC23021B-C844-484C-8E50-ACBAD7E7EE4A.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/E1AAFA19-5E44-8F47-B277-AFB72370E389.root",
    "/eos/home-b/bbapi/CMSSW_13_3_1_patch1/src/NanoAODProduction/NanoAODv13/test/M60-RunIISummer20UL18MiniAODv2/F2510EFC-7848-9944-9EEC-4DEF57B3460F.root"
  ]
}


def process_and_save(dataset_name, files):
    try:
        runner = Runner(
            executor=FuturesExecutor(compression=None, workers=20),
            schema=NanoAODSchema,
            # schemaclass=NanoAODSchema,
            # schemaargs={"version": "v12"},
            savemetrics=True,
        )

        output, _ = runner(
            fileset={dataset_name: files},
            treename="Events",
            processor_instance=HiggsAnalysisProcessor(),
        )

        events = {
            key: ak.Array(val.value) if hasattr(val, "value") else val
            for key, val in output[dataset_name].items()
        }

        if not events or any(len(v) == 0 for v in events.values()):
            print(f"Skipping {dataset_name}: No events")
            return None

        num_events = len(next(iter(events.values())))
        events["dataset"] = ak.Array([dataset_name] * num_events)

        os.makedirs("parquet_files_WH_reco", exist_ok=True)
        out_file = os.path.join("parquet_files_WH_reco", f"{dataset_name}.parquet")
        table = ak.to_arrow_table(events)
        pq.write_table(table, out_file, compression=None)

        return out_file

    except Exception as e:
        print(f"Error processing {dataset_name}: {e}")
        return None


def is_parquet_valid(file_path, delete_if_invalid=False):
    """Check if parquet file exists, is non-empty, and can be opened."""
    try:
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            pq.read_table(file_path)
            return True
    except Exception:
        pass

    if delete_if_invalid and os.path.exists(file_path):
        print(f"Removing broken file: {file_path}")
        try:
            os.remove(file_path)
        except Exception as e:
            print(f"Could not delete {file_path}: {e}")

    return False

fileset_sources = [fileset]

def run_all_datasets(fileset_sources, max_workers=10, max_retries=3):
    failed_datasets = list(fileset_sources[0].keys())  # assume same keys across filesets
    successful = []
    skipped = []

    for attempt in range(1, max_retries + 1):
        if not failed_datasets:
            break

        print(f"\nAttempt {attempt} with {len(failed_datasets)} datasets...")

        # Choose which redirector to use on this attempt
        fileset = fileset_sources[(attempt - 1) % len(fileset_sources)]
        print(f"   → Using redirector set: {['global','fnal','cnaf'][(attempt - 1) % len(fileset_sources)]}")

        remaining = []
        for ds in failed_datasets:
            out_path = os.path.join("parquet_files_WH_reco", f"{ds}.parquet")
            if is_parquet_valid(out_path, delete_if_invalid=True):
                print(f"Skipping {ds} (already processed)")
                skipped.append(ds)
            else:
                remaining.append(ds)

        if not remaining:
            break

        new_failures = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_and_save, name, fileset[name]): name
                for name in remaining if name in fileset
            }

            for future in as_completed(futures):
                dataset = futures[future]
                try:
                    out_file = future.result()
                    if out_file and is_parquet_valid(out_file, delete_if_invalid=True):
                        print(f"Finished {dataset}")
                        successful.append(dataset)
                    else:
                        print(f"No valid file for {dataset}")
                        new_failures.append(dataset)
                except Exception as e:
                    print(f"Error processing {dataset}: {e}")
                    new_failures.append(dataset)

        failed_datasets = new_failures

    if failed_datasets:
        print(f"\nStill failed after retries: {failed_datasets}")

    # Final Summary
    print("\nFinal Summary:")
    print(f"Successful: {len(successful)} → {successful}")
    print(f"Skipped (already valid): {len(skipped)} → {skipped}")
    print(f"Failed after all retries: {len(failed_datasets)} → {failed_datasets}")

    return successful, skipped, failed_datasets


def load_with_label(file):
    try:
        label = os.path.basename(file).split(".")[0]
        arr = ak.from_parquet(file)
        return ak.with_field(arr, label, "dataset")
    except Exception as e:
        print(f"Could not load {file}: {e}")
        return None


# Main execution
# successful, skipped, failed = run_all_datasets(fileset, max_workers=10, max_retries=2)

successful, skipped, failed = run_all_datasets(fileset_sources, max_workers=10, max_retries=3)

# Merge only valid parquet files
parquet_dir = "parquet_files_WH_reco"
files = [f for f in glob.glob(os.path.join(parquet_dir, "*.parquet")) if is_parquet_valid(f)]

if files:
    print(f"\nMerging {len(files)} parquet files...")
    labeled_arrays = [arr for arr in (load_with_label(f) for f in files) if arr is not None]

    if labeled_arrays:
        Events = ak.concatenate(labeled_arrays, axis=0)
        table = ak.to_arrow_table(Events)
        pq.write_table(table, "WH_Reco.parquet", compression=None)
        print("Done! All datasets saved and merged.")
    else:
        print("No valid parquet files to merge.")
else:
    print("No parquet files found for merging.")
