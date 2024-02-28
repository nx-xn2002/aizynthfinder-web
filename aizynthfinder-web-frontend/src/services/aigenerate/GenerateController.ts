// @ts-ignore
/* eslint-disable */
import {request} from '@umijs/max';

/** getImgBySmiles GET /generate/getImgBySmiles */
export async function getImgBySmiles(
  params: API.getImgBySmilesUsingGETParams,
  options?: { [key: string]: any },
) {
  return request<API.BaseResponse>('/generate/getImgBySmiles', {
    method: 'GET',
    params: {
      ...params,
    },
    ...(options || {}),
  });
}

/** getRoutesBySmiles GET /generate/getRoutesBySmiles */
export async function getRoutesBySmiles(
  params: API.getRoutesBySmilesGETParams,
  options?: { [key: string]: any },
) {
  return request<API.BaseResponse>('/generate/getRoutesBySmiles', {
    method: 'GET',
    params: {
      ...params,
    },
    ...(options || {}),
  });
}
