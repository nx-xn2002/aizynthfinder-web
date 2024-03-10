declare namespace API {
  type BaseResponse = {
    code?: number;
    description?: string;
    data;
  };
  type getImgBySmilesUsingGETParams = {
    /** smiles */
    smiles?: string;
  };
  type getRoutesBySmilesGETParams = {
    /** smiles */
    smiles?: string;
    scorers?: string[];
  };
  type getMolsBySmilesGETParams = {
    /** smiles */
    smiles?: string;
  };
  type LoginParams = {
    username?: string;
    password?: string;
    autoLogin?: boolean;
    type?: string;
  };
  type UserInfo = {
    username?: string;
    avatar?: string;
  }
}
