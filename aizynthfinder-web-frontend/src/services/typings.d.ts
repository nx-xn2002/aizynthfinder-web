declare namespace API {
  type BaseResponseString_ = {
    code?: number;
    description?: string;
    data?: string;
  };
  type getImgBySmilesUsingGETParams = {
    /** smiles */
    smiles?: string;
  };
}
