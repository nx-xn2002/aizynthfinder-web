// @ts-ignore
/* eslint-disable */
import {request} from '@umijs/max';

/** getImgBySmiles GET /generate/getImgBySmiles */
export async function getUserByIdUsingGet(
  params: API.getImgBySmilesUsingGETParams,
  options?: { [key: string]: any },
) {
  return request<API.BaseResponseString_>('/generate/getImgBySmiles', {
    method: 'GET',
    params: {
      ...params,
    },
    ...(options || {}),
  });
}
