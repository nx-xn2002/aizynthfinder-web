// @ts-ignore
/* eslint-disable */
import {request} from '@umijs/max';

/** getImgBySmiles GET /generate/getImgBySmiles */
export async function getUserByIdUsingGet(
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

/** 登录接口 POST /user/login */
export async function login(body: API.LoginParams, options?: { [key: string]: any }) {
  return request<API.BaseResponse>('/user/login', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    data: body,
    ...(options || {}),
  });
}

/** 获取当前的用户 GET /user/current */
export async function current(options?: { [key: string]: any }) {
  return request<API.BaseResponse>('/user/current', {
    method: 'GET',
    ...(options || {}),
  });
}

/** 退出登录接口 POST /user/logout */
export async function logout(options?: { [key: string]: any }) {
  return request<Record<string, any>>('/user/logout', {
    method: 'POST',
    ...(options || {}),
  });
}
