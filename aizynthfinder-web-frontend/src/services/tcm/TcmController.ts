// @ts-ignore
/* eslint-disable */
import {request} from '@umijs/max';



/** 查询接口 POST /tcm/search */
export async function listTcmByPageUsingPost(body: API.TcmQueryRequest, options?: { [key: string]: any }) {
  return request<API.BaseResponse>('/tcm/search', {
    method: 'POST',
    withCredentials: true,
    headers: {
      'Content-Type': 'application/json',
    },
    data: body,
    ...(options || {}),
  });
}

/** 查询接口 POST /ingredient/search */
export async function IngredientQueryRequest(body: API.TcmResponse, options?: { [key: string]: any }) {
  return request<API.BaseResponse>('/ingredient/search', {
    method: 'POST',
    withCredentials: true,
    headers: {
      'Content-Type': 'application/json',
    },
    data: body,
    ...(options || {}),
  });
}


